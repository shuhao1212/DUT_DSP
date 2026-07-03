/**
 * user_subband.h —— 8子带多相DFT滤波器组（头文件）
 *
 * 所属项目：项目制实验1 —— 子频带信号分析与合成系统
 * 硬件平台：TMS320C6748 DSP
 *
 * 架构说明：
 *   - 子带数目 D = 8，原型滤波器阶数 N = 64，多相分支长度 L = N/D = 8
 *   - 原型滤波器：Kaiser 窗设计（β = 7.86），运行时由 design_prototype() 生成
 *   - 调制方式：8 点 DFT / IDFT，旋转因子预计算为查找表（LUT）
 *   - 计算效率：相比"调制法"每样点 512 次乘法，多相 DFT 仅需 32 次（提升约 16 倍）
 *
 * 使用步骤：
 *   1. Subband_InitFilterCoeffs()      —— 设计原型滤波器 + 多相分解
 *   2. Subband_InitAnalysisState(&st)  —— 清空分析滤波器组延迟线
 *   3. Subband_InitSynthesisState(&st) —— 清空合成滤波器组延迟线
 *   4. Subband_Pipeline(&ana, &syn, in, out, len) —— 逐帧处理
 *
 * Architecture:
 *   - D = 8 sub-bands, prototype order N = 64, polyphase branch length L = 8
 *   - Kaiser window (beta = 7.86), designed at runtime
 *   - 8-point DFT/IDFT via precomputed LUT
 *
 * Usage:
 *   1. Subband_InitFilterCoeffs()     -- design prototype, decompose polyphase
 *   2. Subband_InitAnalysisState(&st) -- clear delay lines
 *   3. Subband_InitSynthesisState(&st)
 *   4. Subband_Pipeline(&ana, &syn, in, out, len) -- process one frame
 */

#ifndef _USER_SUBBAND_H_
#define _USER_SUBBAND_H_

#include <stdint.h>

/* ======================== 系统常量 ======================== */
#define SUBBAND_D          8                   /* 子带数目 D（同时也是 DFT 点数）         */
#define SUBBAND_ORDER      64                  /* 原型滤波器阶数 N（N = q×D, 取 q=8）    */
#define SUBBAND_L          (SUBBAND_ORDER / SUBBAND_D)  /* 多相分支长度 L = N/D = 8    */
#define SUBBAND_MAX_BLKS   256                 /* 管道缓冲器最大块数（2048/8 = 256）      */

/* ======================== 系统延迟（用于输入输出对齐） ========================
 *  分析→合成滤波器组的群延迟，由原型滤波器阶数 N=64 和多相 DFT 结构决定。
 *  理论值: N-1 = 63 样点
 *  MATLAB 实测（脉冲响应峰值法）: 56 样点
 *
 *  校准方法：触发脉冲测试（设置 pulse_test_trigger = 1 或按 KEY5），
 *  观察 measured_delay 变量得到精确值，然后更新此宏。
 */
#define SUBBAND_SYS_DELAY  56                 /* 系统群延迟（样点数），可校准             */

/* ======================== 脉冲测试（延迟校准用） ======================== */
#define PULSE_POS          100                /* 脉冲注入位置（样点索引）                 */
extern volatile int pulse_test_trigger;       /* 设为 1 触发一次脉冲测试（CCS 或 KEY5）  */
extern volatile int measured_delay;           /* 脉冲测试测得的实际延迟（样点数）         */

/* ======================== 实时性能指标（CCS Expressions 可观测） ======================== */
extern volatile float  error_rms;              /* 每帧 RMS 误差（增益补偿后），越小越好    */
extern volatile float  frame_snr_db;           /* 每帧信噪比 SNR (dB)，越大越好            */
extern volatile int16_t error_max;             /* 每帧峰值误差（绝对值），反映最差情况     */
extern volatile float  subband_energy[8];      /* 8个子带各自的能量，反映频谱分布          */

/* ======================== 误差观测（CCS Graph 调试用） ======================== */
extern int16_t error_buf[2048];               /* 对齐后的输入-输出误差，CCS Graph 可观测  */

/* ======================== 滤波器状态结构体 ======================== */

/**
 * 分析/合成滤波器组状态
 * 
 * 每条分支 (共 D=8 条) 各自维护一个长度为 L=8 的环形延迟缓冲区：
 *   delay[k][0..L-1] —— 第 k 条分支的历史采样值
 *   head[k]          —— 环形缓冲区的写指针（下一次写入位置）
 * 
 * 环形缓冲使用位掩码取模（& (L-1)），因为 L=8 是 2 的幂，避免除法指令。
 */
typedef struct {
    float delay[SUBBAND_D][SUBBAND_L];         /* 环形延迟缓冲区：8分支 × 8阶         */
    int   head[SUBBAND_D];                     /* 写指针（环形，位掩码取模）            */
} Subband_AnaState;

typedef struct {
    float delay[SUBBAND_D][SUBBAND_L];         /* 环形延迟缓冲区：8分支 × 8阶         */
    int   head[SUBBAND_D];                     /* 写指针（环形，位掩码取模）            */
} Subband_SynState;

/* ======================== 公开接口函数 ======================== */

/* 初始化多相滤波器系数（设计 Kaiser 窗原型 + 多相分解，运行时计算） */
void Subband_InitFilterCoeffs(void);

/* 清空分析/合成滤波器组的延迟线状态（上电或复位后调用） */
void Subband_InitAnalysisState(Subband_AnaState *st);
void Subband_InitSynthesisState(Subband_SynState *st);

/* 帧级流水线：分析 → 子带域处理（当前为透传） → 合成
 *   ana, syn : 滤波器状态指针（保持跨帧连续性）
 *   input    : 输入音频帧（int16, 通常 2048 点）
 *   output   : 输出音频帧（int16, 与 input 等长）
 *   len      : 帧长度（必须为 D=8 的整数倍）                       */
void Subband_Pipeline(Subband_AnaState *ana, Subband_SynState *syn,
                      const int16_t input[], int16_t output[], int len);

/* 主示例函数：ADC 采集 → 子带滤波器组 → DAC 输出
 *   按键: KEY1=ADC启, KEY2=ADC停, KEY3=DAC启, KEY4=复位, KEY5=脉冲测试 */
void Subband_Example(void);

#endif /* _USER_SUBBAND_H_ */