/**
 * user_subband.c —— 8子带多相DFT滤波器组（实现文件）
 *
 * 所属项目：项目制实验1 —— 子频带信号分析与合成系统
 * 硬件平台：TMS320C6748 DSP
 *
 * 算法架构：
 *   1. 多相分解：64阶原型滤波器 h[n] 按 D=8 分解为 8条分支 × 8阶 = ana_poly[8][8]
 *   2. 分析块：  8点输入 → 分支卷积 → 8点DFT → 8个复数子带值
 *   3. 合成块：  8个复数子带值 → 8点IDFT → 分支卷积 → 8点输出
 *   4. 帧级流水：ADC帧(2048点) → 256次分析块 → pipe_buf → 256次合成块 → DAC帧
 *
 * 数据流：ADC(乒乓缓冲) → adc_buf → Subband_Pipeline → dac_buf → DAC(乒乓缓冲)
 *
 * Architecture: standard polyphase DFT, N=64 prototype, L=8, Kaiser window
 * Uses ADC/DAC dual-flag + ping-pong buffer, compatible with existing project
 */

#include "adc_api.h"
#include "dac_api.h"
#include "system.h"
#include "delay.h"
#include "led_api.h"
#include "key_api.h"
#include "user_subband.h"
#include "math.h"
#include "string.h"
#include <stdint.h>

#define PI 3.14159265358979323846f

/* ========================================================================
 *  误差观测缓冲（全局，CCS Graph 工具可实时显示）
 *
 *  error_buf[i] = 对齐后的 输入 - 输出 之差
 *  由于系统存在 SUBBAND_SYS_DELAY 个样点的群延迟，直接相减无效。
 *  对齐方式：error_buf[i] = adc_buf[i - DELAY] - dac_buf[i]  (i >= DELAY)
 *  前 DELAY 个样点置零（无对应输入，属于滤波器暂态响应区间）。
 *
 *  在 CCS 中添加到 Graph：选 error_buf, 长度 2048, 类型 16-bit signed int
 * ======================================================================== */
int16_t error_buf[2048];

/* ========================================================================
 *  脉冲测试变量（延迟校准用，CCS Watch 窗口可观测）
 *
 *  使用方法：
 *    方法1（CCS）: 在 CCS Expressions 窗口将 pulse_test_trigger 设为 1
 *    方法2（按键）: 按 KEY5（此时 LED 切换暂停，优先执行脉冲测试）
 *
 *  执行后 measured_delay 显示测得的系统延迟，将其更新到
 *  user_subband.h 的 SUBBAND_SYS_DELAY 宏中即可。
 *  同时 dac_buf 和 error_buf 会输出脉冲响应波形，可在示波器或 Graph 上观察。
 * ======================================================================== */
volatile int pulse_test_trigger = 0;            /* 触发标志：设为 1 执行一次            */
volatile int measured_delay = 0;                /* 测量结果：脉冲输出峰值位置 - PULSE_POS */

/* ========================================================================
 *  实时性能指标（全局 volatile，CCS Expressions 窗口实时刷新）
 *  每帧更新一次，可直接拖入 CCS Watch 窗口观察数值变化。
 * ======================================================================== */
volatile float  error_rms = 0.0f;               /* 每帧 RMS 误差（增益补偿后）            */
volatile float  frame_snr_db = 0.0f;            /* 每帧信噪比 (dB)                        */
volatile int16_t error_max = 0;                 /* 每帧峰值误差（绝对值）                 */
volatile float  subband_energy[8] = {0};        /* 8个子带能量分布                        */

/* ========================================================================
 *  8点 DFT 旋转因子查找表（预计算，避免运行时调用 sin/cos）
 *
 *  W8_COS[k][n] = cos(2π·k·n / D)    —— DFT 旋转因子实部
 *  W8_SIN[k][n] = sin(2π·k·n / D)    —— DFT 旋转因子虚部
 *
 *  其中 k = 0..7（频域索引），n = 0..7（时域索引），D = 8
 *
 *  设计考量：
 *    - 8×8×2 = 128 个 float，仅占 512 字节 Flash
 *    - 免去实时 cosf/sinf 调用，C6748 单周期即可完成查表+乘法
 *    - DFT 公式: X[k] = Σ x[n]·exp(-j·2π·k·n/D)
 * ======================================================================== */
/* cos(2*PI*k*n/8), sin(2*PI*k*n/8) for k,n = 0..7 */
static const float W8_COS[SUBBAND_D][SUBBAND_D] = {
    { 1.000000f, 1.000000f, 1.000000f, 1.000000f, 1.000000f, 1.000000f, 1.000000f, 1.000000f },
    { 1.000000f, 0.707107f, 0.000000f,-0.707107f,-1.000000f,-0.707107f,-0.000000f, 0.707107f },
    { 1.000000f, 0.000000f,-1.000000f, 0.000000f, 1.000000f, 0.000000f,-1.000000f, 0.000000f },
    { 1.000000f,-0.707107f,-0.000000f, 0.707107f,-1.000000f, 0.707107f, 0.000000f,-0.707107f },
    { 1.000000f,-1.000000f, 1.000000f,-1.000000f, 1.000000f,-1.000000f, 1.000000f,-1.000000f },
    { 1.000000f,-0.707107f, 0.000000f, 0.707107f,-1.000000f, 0.707107f,-0.000000f,-0.707107f },
    { 1.000000f, 0.000000f,-1.000000f, 0.000000f, 1.000000f, 0.000000f,-1.000000f, 0.000000f },
    { 1.000000f, 0.707107f, 0.000000f,-0.707107f,-1.000000f,-0.707107f, 0.000000f, 0.707107f },
};
static const float W8_SIN[SUBBAND_D][SUBBAND_D] = {
    { 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f },
    { 0.000000f, 0.707107f, 1.000000f, 0.707107f, 0.000000f,-0.707107f,-1.000000f,-0.707107f },
    { 0.000000f, 1.000000f, 0.000000f,-1.000000f, 0.000000f, 1.000000f, 0.000000f,-1.000000f },
    { 0.000000f, 0.707107f,-1.000000f, 0.707107f, 0.000000f,-0.707107f, 1.000000f,-0.707107f },
    { 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.000000f },
    { 0.000000f,-0.707107f, 1.000000f,-0.707107f, 0.000000f, 0.707107f,-1.000000f, 0.707107f },
    { 0.000000f,-1.000000f, 0.000000f, 1.000000f, 0.000000f,-1.000000f, 0.000000f, 1.000000f },
    { 0.000000f,-0.707107f,-1.000000f,-0.707107f, 0.000000f, 0.707107f, 1.000000f, 0.707107f },
};

/* ========================================================================
 *  多相滤波器系数（全局静态，由 Subband_InitFilterCoeffs() 初始化）
 *
 *  ana_poly[k][l] = h[k + l·D]          —— 分析滤波器组的多相分量
 *  syn_poly[k][l] = h[k + (L-1-l)·D]    —— 合成滤波器组的多相分量（时间翻转）
 *
 *  原型滤波器 h[n] 通过 Kaiser 窗法设计（β=7.86, fc=1/(2D)=0.0625）
 *  然后在运行时完成多相分解。
 *
 *  注意：syn_poly 是 ana_poly 的时间翻转版本，
 *        这是为了实现"分析+合成=完美重建"所必需的。
 * ======================================================================== */
static float ana_poly[SUBBAND_D][SUBBAND_L];   /* 分析: e[k][l] = h[k + l*D]      */
static float syn_poly[SUBBAND_D][SUBBAND_L];   /* 合成: 时间翻转版本               */

/* ========================================================================
 *  管道缓冲器（帧级中间存储，static 避免栈溢出）
 *
 *  pipe_buf[子带索引 k][块索引 blk][0/1=实部/虚部]
 *  最大容量: 8子带 × 256块 × 2 = 4096 个 float ≈ 16KB
 *
 *  作用：在分析阶段和合成阶段之间缓存子带域数据，
 *        为后续插入子带域处理（如增益控制）预留接口。
 * ======================================================================== */
static float pipe_buf[SUBBAND_D][SUBBAND_MAX_BLKS][2];

/* ========================================================================
 *  第零阶修正贝塞尔函数 I₀(x) —— Kaiser 窗的核心函数
 *
 *  级数展开: I₀(x) = Σ(t=0→∞) [(x/2)^t / t!]²
 *
 *  实现细节：
 *    - 迭代最多 30 项，精度达到 1e-14（单精度浮点已足够）
 *    - 提前退出条件：当前项的平方 < 累加和 × 1e-14
 *    - 对应 MATLAB 中的 besseli(0, x)
 * ======================================================================== */
static float bessel_i0(float x) {
    float t = 1.0f, s = 1.0f, xh = x * 0.5f;  /* t=第n项, s=累加和, xh=x/2 */
    int i;
    for (i = 1; i < 30; i++) {
        t *= xh / (float)i;                    /* 递推计算 (x/2)^i / i!      */
        float a = t * t;                       /* 当前项的平方               */
        s += a;                                /* 累加                       */
        if (a < s * 1.0e-14f) break;           /* 收敛判定，提前退出         */
    }
    return s;
}

/* ========================================================================
 *  Kaiser 窗原型低通滤波器设计（运行时计算，非硬编码）
 *
 *  参数说明：
 *    h[]   —— 输出，长度为 N 的滤波器系数数组
 *    N     —— 滤波器阶数（本项目固定为 64）
 *    fc    —— 归一化截止频率（本项目为 1/(2D) = 0.0625，即 π/16 rad）
 *    beta  —— Kaiser 窗形状参数（β=7.86 对应阻带衰减 ≈ 60dB）
 *
 *  设计公式：
 *    h_ideal[n] = sin(ωc·(n-α)) / (π·(n-α))     —— 理想低通 sinc 函数
 *    w[n]       = I₀(β·√(1-(2n/(N-1)-1)²)) / I₀(β)  —— Kaiser 窗
 *    h[n]       = h_ideal[n] · w[n]              —— 加窗
 *    h[n]      /= Σ h[n]                         —— 直流增益归一化
 *
 *  其中 α = (N-1)/2 为对称中心。
 *
 *  设计考量：
 *    - 在 DSP 上电时现场计算，方便改参数
 *    - Kaiser 窗阻带衰减远超 Hanning 窗（82dB vs 44dB）
 *    - N=64 = D×8，保证多相分解时各分支长度一致
 * ======================================================================== */
static void design_prototype(float h[], int N, float fc, float beta) {
    int n;
    float alpha = (N - 1) * 0.5f;               /* 滤波器对称中心（线性相位） */
    float denom = bessel_i0(beta);               /* I₀(β)，窗函数分母           */
    float wc = 2.0f * PI * fc;                   /* 数字截止角频率 ωc           */
    float sum = 0.0f;

    for (n = 0; n < N; n++) {
        float delta = n - alpha;                 /* 到对称中心的距离            */
        float ideal;
        /* sinc 函数: sin(ωc·t) / (π·t)，t=0 时取极限 ωc/π = 2fc */
        if (fabsf(delta) < 1.0e-8f)
            ideal = 2.0f * fc;
        else
            ideal = sinf(wc * delta) / (PI * delta);

        /* Kaiser 窗: I₀(β·√(1-(2n/(N-1)-1)²)) / I₀(β) */
        float arg = (2.0f * n / (N - 1)) - 1.0f;
        arg = 1.0f - arg * arg;
        if (arg < 0.0f) arg = 0.0f;
        float win = bessel_i0(beta * sqrtf(arg)) / denom;

        h[n] = ideal * win;                      /* 理想响应 × 窗函数          */
        sum += h[n];
    }
    /* 直流增益归一化：使 Σ h[n] = 1.0，保证分析/合成级联增益为 1 */
    for (n = 0; n < N; n++) h[n] /= sum;
}

/* ========================================================================
 *  多相分解 —— 将原型滤波器 h[n] 分解为 D 条长度为 L 的子滤波器
 *
 *  分解公式：
 *    分析多相分量: ana_poly[k][l] = h[k + l·D]      (k=分支, l=抽头)
 *    合成多相分量: syn_poly[k][l] = h[k + (L-1-l)·D]  (时间翻转)
 *
 *  合成滤波器的时间翻转是为了实现"分析+合成 = 近似完美重建"。
 *
 *  注意：此函数必须在系统初始化时调用一次，且必须在
 *        Subband_InitAnalysisState / Subband_InitSynthesisState 之前。
 * ======================================================================== */
void Subband_InitFilterCoeffs(void) {
    float h[SUBBAND_ORDER];                      /* 原型滤波器系数（临时）      */
    int k, l;

    /* 第1步：Kaiser 窗设计 64 阶原型低通滤波器
     *   fc = 1/(2D) = 0.0625（归一化频率），β = 7.86（阻带衰减 ≈ 60dB） */
    design_prototype(h, SUBBAND_ORDER, 1.0f / (2.0f * SUBBAND_D), 7.86f);

    /* 第2步：多相分解
     *   ana_poly[k][l] = h[k + l*D]          —— 分析滤波器多相分量
     *   syn_poly[k][l] = h[k + (L-1-l)*D]    —— 合成滤波器多相分量（翻转） */
    for (k = 0; k < SUBBAND_D; k++) {
        for (l = 0; l < SUBBAND_L; l++) {
            ana_poly[k][l] = h[k + l * SUBBAND_D];
            syn_poly[k][l] = h[k + (SUBBAND_L-1-l) * SUBBAND_D];
        }
    }
}

/* ========================================================================
 *  滤波器状态初始化 —— 清空所有延迟线
 *
 *  将每条分支的环形缓冲区全部置零，写指针复位到 0。
 *  在以下场景调用：
 *    - 系统上电初始化
 *    - KEY4 按下（用户手动复位滤波器）
 *    - 处理不同音频段落之间（避免上一段尾部影响下一段头部）
 * ======================================================================== */
void Subband_InitAnalysisState(Subband_AnaState *st) {
    int k, t;
    for (k = 0; k < SUBBAND_D; k++) {
        for (t = 0; t < SUBBAND_L; t++) st->delay[k][t] = 0.0f;
        st->head[k] = 0;                         /* 写指针归零                  */
    }
}
void Subband_InitSynthesisState(Subband_SynState *st) {
    int k, t;
    for (k = 0; k < SUBBAND_D; k++) {
        for (t = 0; t < SUBBAND_L; t++) st->delay[k][t] = 0.0f;
        st->head[k] = 0;
    }
}

/* ========================================================================
 *  8 点 DFT（离散傅里叶变换） —— 基于预计算查找表
 *
 *  公式：X[k] = Σ(n=0..7) x[n] · exp(-j·2π·k·n/8)
 *  展开：re_out[k] = Σ( re_in[n]·cos + im_in[n]·sin )
 *        im_out[k] = Σ( im_in[n]·cos - re_in[n]·sin )
 *
 *  计算量：8×8 = 64 次乘法 + 64 次加法（每 8 个输入样点调用一次，
 *          平摊到每样点仅 8 次乘法，效率远高于直接 DFT 的 O(D²)）
 * ======================================================================== */
static void dft8(const float re_in[], const float im_in[], float re_out[], float im_out[]) {
    int k, n;
    for (k = 0; k < SUBBAND_D; k++) {            /* 遍历频域索引 k=0..7        */
        float sr = 0, si = 0;
        for (n = 0; n < SUBBAND_D; n++) {        /* 遍历时域索引 n=0..7        */
            float c = W8_COS[k][n], s = W8_SIN[k][n];
            /* 复数乘法: (re+j·im) × (cos-j·sin) = (re·cos+im·sin) + j(im·cos-re·sin) */
            sr += re_in[n] * c + im_in[n] * s;
            si += im_in[n] * c - re_in[n] * s;
        }
        re_out[k] = sr; im_out[k] = si;
    }
}

/* ========================================================================
 *  8 点 IDFT（离散傅里叶逆变换） —— 基于预计算查找表
 *
 *  公式：x[n] = (1/D) · Σ(k=0..7) X[k] · exp(+j·2π·k·n/8)
 *  展开：re_out[n] = (1/D)·Σ( re_in[k]·cos - im_in[k]·sin )
 *        im_out[n] = (1/D)·Σ( re_in[k]·sin + im_in[k]·cos )
 *
 *  注意：除以 D 是 IDFT 的固有缩放，DFT→IDFT 级联恢复原始信号。
 *        在分析→合成链中，DFT 不缩放，IDFT 除以 D，这是标准做法。
 * ======================================================================== */
static void idft8(const float re_in[], const float im_in[], float re_out[], float im_out[]) {
    int k, n;
    for (n = 0; n < SUBBAND_D; n++) {            /* 遍历时域索引 n=0..7        */
        float sr = 0, si = 0;
        for (k = 0; k < SUBBAND_D; k++) {        /* 遍历频域索引 k=0..7        */
            float c = W8_COS[k][n], s = W8_SIN[k][n];
            /* 复数乘法: (re+j·im) × (cos+j·sin) = (re·cos-im·sin) + j(re·sin+im·cos) */
            sr += re_in[k] * c - im_in[k] * s;
            si += re_in[k] * s + im_in[k] * c;
        }
        re_out[n] = sr / SUBBAND_D;              /* IDFT 缩放因子 1/D          */
        im_out[n] = si / SUBBAND_D;
    }
}

/* ========================================================================
 *  分析块 —— 将 8 个时域采样点分解为 8 个复数子带值
 *
 *  参数：
 *    st  —— 分析滤波器状态（含延迟线，跨块保持）
 *    in  —— 输入，8 个 int16 时域采样点
 *    out —— 输出，8 个复数子带值 out[k][0]=实部, out[k][1]=虚部
 *
 *  处理流程（与多相 DFT 理论严格对应）：
 *    第1步：8 点输入分别写入 8 条分支的环形延迟线
 *    第2步：每条分支执行 L=8 阶 FIR 卷积（多相滤波）
 *    第3步：8 路卷积结果送入 8 点 DFT，输出即为 8 路子带复数信号
 *
 *  环形缓冲区实现细节：
 *    - head[k] 始终指向"下一次写入位置"
 *    - tap=0 读取最新的采样值（head-1），tap=L-1 读取最旧的（head）
 *    - 取模用位掩码 & (L-1)，因 L=8 是 2 的幂，避免除法
 * ======================================================================== */
void Subband_AnalysisBlock(Subband_AnaState *st, const int16_t in[SUBBAND_D], float out[SUBBAND_D][2]) {
    int k, tap;
    float u_re[SUBBAND_D], u_im[SUBBAND_D];      /* 多相卷积输出（8路）        */
    float dft_re[SUBBAND_D], dft_im[SUBBAND_D];  /* DFT 输出（8路复数）        */

    for (k = 0; k < SUBBAND_D; k++) {
        /* --- 写入输入采样到环形缓冲区 --- */
        st->delay[k][st->head[k]] = (float)in[k];
        st->head[k] = (st->head[k] + 1) & (SUBBAND_L - 1);  /* 位掩码取模 */

        /* --- L=8 阶多相 FIR 卷积 --- */
        float acc = 0;
        for (tap = 0; tap < SUBBAND_L; tap++) {
            /* 从最新到最旧依次读取: head-1, head-2, ..., head-L */
            int idx = (st->head[k] + SUBBAND_L - 1 - tap) & (SUBBAND_L - 1);
            acc += ana_poly[k][tap] * st->delay[k][idx];
        }
        u_re[k] = acc;                             /* 分支 k 的卷积结果（实数） */
        u_im[k] = 0;                               /* 输入为实数，虚部为 0      */
    }

    /* --- 8 点 DFT：将多相域变换到子带域 --- */
    dft8(u_re, u_im, dft_re, dft_im);

    /* --- 输出：8 个复数子带值 --- */
    for (k = 0; k < SUBBAND_D; k++) {
        out[k][0] = dft_re[k]; out[k][1] = dft_im[k];
    }
}

/* ========================================================================
 *  合成块 —— 将 8 个复数子带值重建为 8 个时域采样点
 *
 *  参数：
 *    st  —— 合成滤波器状态（含延迟线）
 *    in  —— 输入，8 个复数子带值 in[k][0]=实部, in[k][1]=虚部
 *    out —— 输出，8 个 int16 重建时域采样点
 *
 *  处理流程（分析块的逆过程）：
 *    第1步：8 个复数子带值送入 8 点 IDFT
 *    第2步：IDFT 的 8 路输出分别进入 8 条分支的延迟线
 *    第3步：每条分支执行 L=8 阶 FIR 卷积（合成多相滤波）
 *    第4步：卷积结果饱和钳位到 int16 范围后输出
 *
 *  饱和处理：防止浮点累积误差导致的超出 int16 范围，
 *            避免 DAC 输出出现"翻转"（32767 → -32768）的刺耳噪声。
 * ======================================================================== */
void Subband_SynthesisBlock(Subband_SynState *st, const float in[SUBBAND_D][2], int16_t out[SUBBAND_D]) {
    int k, tap;
    float w_re[SUBBAND_D], w_im[SUBBAND_D];      /* IDFT 输入（8路复数）       */
    float idft_re[SUBBAND_D], idft_im[SUBBAND_D];/* IDFT 输出（8路）           */

    /* --- 准备 IDFT 输入 --- */
    for (k = 0; k < SUBBAND_D; k++) {
        w_re[k] = in[k][0]; w_im[k] = in[k][1];
    }

    /* --- 8 点 IDFT：将子带域变换回多相域 --- */
    idft8(w_re, w_im, idft_re, idft_im);

    for (k = 0; k < SUBBAND_D; k++) {
        /* --- 写入 IDFT 输出到环形缓冲区 --- */
        st->delay[k][st->head[k]] = idft_re[k];    /* 只需实部（虚部为 0）     */
        st->head[k] = (st->head[k] + 1) & (SUBBAND_L - 1);

        /* --- L=8 阶多相 FIR 卷积 --- */
        float acc = 0;
        for (tap = 0; tap < SUBBAND_L; tap++) {
            int idx = (st->head[k] + SUBBAND_L - 1 - tap) & (SUBBAND_L - 1);
            acc += syn_poly[k][tap] * st->delay[k][idx];
        }

        /* --- 饱和钳位到 int16 范围 [-32768, 32767] --- */
        if (acc >  32767.0f) acc =  32767.0f;
        if (acc < -32768.0f) acc = -32768.0f;
        out[k] = (int16_t)acc;                     /* 分支 k 的重建输出         */
    }
}

/* ========================================================================
 *  帧级流水线 —— 处理一帧完整的音频数据
 *
 *  流程：分析所有块 → [子带域处理（当前为透传）] → 合成所有块
 *
 *  参数：
 *    ana, syn —— 滤波器状态（跨帧保持连续性）
 *    input    —— 输入音频帧（int16 数组，通常 2048 点）
 *    output   —— 输出音频帧（int16 数组，与 input 等长）
 *    len      —— 帧长度（必须为 D=8 的整数倍）
 *
 *  两阶段设计的意义：
 *    正向路径（分析）：逐块分析，将所有子带值存入 pipe_buf
 *    反向路径（合成）：从 pipe_buf 读取子带值，逐块合成
 *    
 *    这种设计为"子带域处理"预留了接口——在两次循环之间对
 *    pipe_buf 中的子带数据做任意处理（如每子带增益控制），
 *    即可实现图示均衡器等扩展功能。
 *
 *  时序约束：C6748 @ 456MHz，2048 点 Pipeline 应在 40.96ms 帧周期内完成。
 * ======================================================================== */
void Subband_Pipeline(Subband_AnaState *ana, Subband_SynState *syn,
                      const int16_t input[], int16_t output[], int len) {
    int blk_cnt = len / SUBBAND_D;               /* 块数 = 2048/8 = 256         */
    int k, blk;
    float sbuf[SUBBAND_D][2];                    /* 临时子带缓冲（单个块）      */

    if (blk_cnt > SUBBAND_MAX_BLKS) blk_cnt = SUBBAND_MAX_BLKS;

    /* ===== 正向路径：分析滤波器组 ===== */
    /* 逐块分析 256 次，将 8×256 个子带值存入管道缓冲 */
    for (blk = 0; blk < blk_cnt; blk++) {
        Subband_AnalysisBlock(ana, &input[blk * SUBBAND_D], sbuf);
        for (k = 0; k < SUBBAND_D; k++) {
            pipe_buf[k][blk][0] = sbuf[k][0];    /* 子带 k, 块 blk, 实部      */
            pipe_buf[k][blk][1] = sbuf[k][1];    /* 子带 k, 块 blk, 虚部      */
        }
    }

    /* ===== 子带域处理（扩展点） ===== */
    /* 当前为"透传"模式，不修改 pipe_buf。
     * 若需实现图示均衡器，在此处对 pipe_buf[k][blk][*] 乘以增益系数即可。 */

    /* ===== 子带能量统计（诊断用，存入 subband_energy[]） ===== */
    {
        int k;
        for (k = 0; k < SUBBAND_D; k++) {
            float energy = 0.0f;
            int blk;
            for (blk = 0; blk < blk_cnt; blk++) {
                float re = pipe_buf[k][blk][0];
                float im = pipe_buf[k][blk][1];
                energy += re * re + im * im;
            }
            subband_energy[k] = energy / (float)blk_cnt;  /* 每块平均能量 */
        }
    }

    /* ===== 反向路径：合成滤波器组 ===== */
    /* 逐块从管道缓冲读取子带值，合成 256 块输出 */
    for (blk = 0; blk < blk_cnt; blk++) {
        for (k = 0; k < SUBBAND_D; k++) {
            sbuf[k][0] = pipe_buf[k][blk][0];
            sbuf[k][1] = pipe_buf[k][blk][1];
        }
        Subband_SynthesisBlock(syn, sbuf, &output[blk * SUBBAND_D]);
    }
}

/* ========================================================================
 *  主示例函数 —— ADC 采集 → 子带滤波器组 → DAC 输出
 *
 *  完整数据流：
 *    ADC (50kHz, 乒乓缓冲) → adc_buf[2048] → Subband_Pipeline → dac_buf[2048] → DAC
 *
 *  按键交互功能：
 *    KEY1 —— 启动 ADC 采集
 *    KEY2 —— 停止 ADC 采集
 *    KEY3 —— 启动 DAC 输出
 *    KEY4 —— 复位滤波器状态（停止 ADC/DAC → 清空延迟线 → 重启）
 *    KEY5 —— 触发脉冲测试（延迟校准，红灯亮表示完成）
 *
 *  蓝色 LED 在每次成功完成一帧处理时翻转，用于观察帧处理节奏。
 *
 *  初始化顺序说明（严格依赖）：
 *    系统 → LED → 按键 → 滤波器系数 → 滤波器状态 → ADC → DAC → 启动
 *
 *  乒乓缓冲机制：
 *    ADC 和 DAC 各有两个缓冲区（PONG/PING），硬件在采集/输出一个缓冲区
 *    的同时，CPU 可以处理另一个缓冲区的数据，实现零拷贝流水线。
 * ======================================================================== */
void Subband_Example(void) {
    unsigned char ad_done = 0;                     /* ADC 数据就绪标志           */
    static Subband_AnaState ana;                   /* 分析滤波器状态（静态）     */
    static Subband_SynState syn;                   /* 合成滤波器状态（静态）     */
    static int16_t adc_buf[2048], dac_buf[2048];  /* ADC/DAC 数据缓冲（静态）   */

    /* ===== 系统初始化 ===== */
    Sys_Init();
    Led_Init();
    Key_Init();
    Led_Control(LED1_CORE, LED_ON);               /* 点亮核心 LED 表示系统运行  */
    Led_Control(LED2_CORE, LED_ON);

    /* ===== 滤波器初始化 ===== */
    Subband_InitFilterCoeffs();                    /* 设计原型 + 多相分解        */
    Subband_InitAnalysisState(&ana);               /* 清空分析延迟线              */
    Subband_InitSynthesisState(&syn);              /* 清空合成延迟线              */

    /* ===== ADC/DAC 初始化 ===== */
    Adc_Init(ADC_50KHZ, 2048);                     /* ADC: 50kHz 采样, 2048点/帧 */
    Dac_Init(DAC_50KHZ, 2048, DAC_CHANNEL_12);    /* DAC: 50kHz, 双通道输出     */
    Adc_Start();                                   /* 启动 ADC 乒乓采集           */
    Dac_Start();                                   /* 启动 DAC 乒乓输出           */

    /* ===== 主循环：乒乓缓冲 + 按键轮询 ===== */
    while (1) {
        /* --- ADC 数据就绪：从乒乓缓冲拷贝到 adc_buf --- */
        if (FLAG_AD == 1) {
            FLAG_AD = 0;
            ad_done = 1;

            if (AD_Ping_Pong == AD_BUFFER_PONG)
                memcpy(adc_buf, AD_CH1_Buf0, sizeof(int16_t) * 2048);
            else
                memcpy(adc_buf, AD_CH1_Buf1, sizeof(int16_t) * 2048);
        }

        /* --- DAC 就绪 + ADC 有数据：执行一帧处理 --- */
        if (FLAG_DA == 1 && ad_done == 1) {
            FLAG_DA = 0;
            ad_done = 0;                           /* 消费 ADC 数据               */

            /* ===== 脉冲测试模式（延迟校准） =====
             *  触发后，用人工脉冲替代真实 ADC 数据，测量系统实际延迟。
             *  测得的延迟存入 measured_delay，可查看后更新 SUBBAND_SYS_DELAY。
             *
             *  脉冲信号：除 PULSE_POS=100 处为 32767 外全为零。
             *  输出脉冲响应会同时出现在 dac_buf（示波器）和 error_buf（Graph）中。
             */
            if (pulse_test_trigger) {
                int i;
                pulse_test_trigger = 0;            /* 清除触发标志（单次执行）    */

                /* 用脉冲替代 adc_buf */
                for (i = 0; i < 2048; i++) adc_buf[i] = 0;
                adc_buf[PULSE_POS] = 32767;        /* 正脉冲，int16 满幅         */

                /* 正常运行 Pipeline */
                Subband_Pipeline(&ana, &syn, adc_buf, dac_buf, 2048);

                /* 在 dac_buf 中搜索输出峰值位置（跳过前几个暂态点） */
                {
                    int peak_idx = PULSE_POS + 10;  /* 从稍后位置开始搜索       */
                    int16_t peak_val = 0;
                    for (i = PULSE_POS + 10; i < 2048; i++) {
                        int16_t abs_val = (dac_buf[i] >= 0) ? dac_buf[i]
                                                            : (int16_t)(-dac_buf[i]);
                        if (abs_val > peak_val) {
                            peak_val = abs_val;
                            peak_idx = i;
                        }
                    }
                    measured_delay = peak_idx - PULSE_POS;
                }

                /* 误差缓冲也填入脉冲响应（方便 Graph 观察脉冲形状） */
                for (i = 0; i < 2048; i++) {
                    error_buf[i] = dac_buf[i];     /* 直接显示脉冲响应波形       */
                }

                Led_Control(LED1_RED, LED_ON);     /* 红灯亮：脉冲测试完成       */
                Led_Control(LED2_RED, LED_ON);
            } else {
                /* ===== 正常模式：处理真实 ADC 数据 ===== */

                /* 核心：子带分析 → 透传 → 子带合成 */
                Subband_Pipeline(&ana, &syn, adc_buf, dac_buf, 2048);

                /* ===== 误差计算（延迟对齐 + 增益补偿） =====
                 *  系统群延迟 = SUBBAND_SYS_DELAY 样点。
                 *  系统增益 ≈ 1/D² = 1/64（多相 DFT 固有衰减）。
                 *
                 *  补偿方式：将输出放大 D²=64 倍后再与延迟后的输入比较。
                 *  使用 int32 中间变量防止溢出。
                 */
                {
                    int i;
                    /* 暂态区清零 */
                    for (i = 0; i < SUBBAND_SYS_DELAY; i++) {
                        error_buf[i] = 0;
                    }

                    /* 有效区：增益补偿后相减 */
                    float sum_sq_err = 0.0f, sum_sq_sig = 0.0f;
                    int16_t peak = 0;
                    int valid_cnt = 2048 - SUBBAND_SYS_DELAY;

                    for (i = SUBBAND_SYS_DELAY; i < 2048; i++) {
                        int32_t err = (int32_t)adc_buf[i - SUBBAND_SYS_DELAY]
                                    - (int32_t)dac_buf[i] * SUBBAND_D * SUBBAND_D;
                        if (err >  32767) err =  32767;
                        if (err < -32768) err = -32768;
                        error_buf[i] = (int16_t)err;

                        /* 累加用于 RMS / SNR / 峰值 计算 */
                        float ef = (float)err;
                        float sf = (float)adc_buf[i - SUBBAND_SYS_DELAY];
                        sum_sq_err += ef * ef;
                        sum_sq_sig += sf * sf;

                        int16_t abs_err = (err >= 0) ? (int16_t)err : (int16_t)(-err);
                        if (abs_err > peak) peak = abs_err;
                    }

                    /* 更新全局性能指标 */
                    if (valid_cnt > 0) {
                        error_rms = sqrtf(sum_sq_err / (float)valid_cnt);
                        if (sum_sq_err > 1.0e-12f && sum_sq_sig > 1.0e-12f) {
                            frame_snr_db = 10.0f * log10f(sum_sq_sig / sum_sq_err);
                        } else {
                            frame_snr_db = 99.0f;     /* 信号太弱，SNR 无意义    */
                        }
                    } else {
                        error_rms = 0.0f;
                        frame_snr_db = 0.0f;
                    }
                    error_max = peak;
                }

                Led_Control(LED1_BLUE, LED_TOGGLE);/* 蓝色LED翻转：帧处理指示     */
            }

            /* 将处理结果拷贝到 DAC 乒乓缓冲 */
            if (DA_Ping_Pong == DA_BUFFER_PONG) {
                memcpy(DA_CH1_Buf0, dac_buf, sizeof(int16_t) * 2048);
                memcpy(DA_CH2_Buf0, dac_buf, sizeof(int16_t) * 2048);
            } else {
                memcpy(DA_CH1_Buf1, dac_buf, sizeof(int16_t) * 2048);
                memcpy(DA_CH2_Buf1, dac_buf, sizeof(int16_t) * 2048);
            }
        }

        /* ===== 按键处理 ===== */

        /* KEY1: 启动 ADC 采集 */
        if (FLAG_KEY1 == 1) { FLAG_KEY1 = 0; Adc_Start(); }

        /* KEY2: 停止 ADC 采集 */
        if (FLAG_KEY2 == 1) { FLAG_KEY2 = 0; Adc_Stop();  }

        /* KEY3: 启动 DAC 输出 */
        if (FLAG_KEY3 == 1) { FLAG_KEY3 = 0; Dac_Start(); }

        /* KEY4: 复位滤波器（停止 → 清空延迟线 → 重启） */
        if (FLAG_KEY4 == 1) {
            FLAG_KEY4 = 0;
            Dac_Stop(); Adc_Stop();                /* 先停 ADC/DAC 避免数据错位  */
            Subband_InitAnalysisState(&ana);       /* 清空所有分析延迟线          */
            Subband_InitSynthesisState(&syn);      /* 清空所有合成延迟线          */
            Adc_Start(); Dac_Start();              /* 重新启动                    */
        }

        /* KEY5: 触发脉冲测试（延迟校准） */
        if (FLAG_KEY5 == 1) {
            FLAG_KEY5 = 0;
            pulse_test_trigger = 1;                /* 下一帧执行脉冲测试          */
        }
    }
}