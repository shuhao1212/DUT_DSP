/**
 * user_subband.c - 8-band Polyphase DFT Filter Bank (adapted)
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

/* ========== Precomputed 8-point DFT Twiddle Factor LUT ========== */
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

/* ========== Polyphase Filter Coefficients ========== */
static float ana_poly[SUBBAND_D][SUBBAND_L];   /* analysis:  e[k][l] = h[k + l*D]       */
static float syn_poly[SUBBAND_D][SUBBAND_L];   /* synthesis: branch time-reversed        */

/* ========== Pipeline Buffer (static, not on stack) ========== */
static float pipe_buf[SUBBAND_D][SUBBAND_MAX_BLKS][2];

/* ========== Bessel I0 (for Kaiser window) ========== */
static float bessel_i0(float x) {
    float t = 1.0f, s = 1.0f, xh = x * 0.5f;
    int i;
    for (i = 1; i < 30; i++) {
        t *= xh / (float)i;
        float a = t * t;
        s += a;
        if (a < s * 1.0e-14f) break;
    }
    return s;
}

/* ========== Kaiser Window Prototype Design (N=64, fc=1/(2D)=0.0625, beta=7.86) ========== */
static void design_prototype(float h[], int N, float fc, float beta) {
    int n;
    float alpha = (N - 1) * 0.5f;
    float denom = bessel_i0(beta);
    float wc = 2.0f * PI * fc;
    float sum = 0.0f;

    for (n = 0; n < N; n++) {
        float delta = n - alpha;
        float ideal;
        if (fabsf(delta) < 1.0e-8f)
            ideal = 2.0f * fc;
        else
            ideal = sinf(wc * delta) / (PI * delta);

        float arg = (2.0f * n / (N - 1)) - 1.0f;
        arg = 1.0f - arg * arg;
        if (arg < 0.0f) arg = 0.0f;
        float win = bessel_i0(beta * sqrtf(arg)) / denom;

        h[n] = ideal * win;
        sum += h[n];
    }
    for (n = 0; n < N; n++) h[n] /= sum;   /* normalize to unity DC gain */
}

/* ========== Polyphase Decomposition ========== */
void Subband_InitFilterCoeffs(void) {
    float h[SUBBAND_ORDER];
    int k, l;

    design_prototype(h, SUBBAND_ORDER, 1.0f / (2.0f * SUBBAND_D), 7.86f);

    for (k = 0; k < SUBBAND_D; k++) {
        for (l = 0; l < SUBBAND_L; l++) {
            ana_poly[k][l] = h[k + l * SUBBAND_D];
            syn_poly[k][l] = h[k + (SUBBAND_L-1-l) * SUBBAND_D];
        }
    }
}

/* ========== State Initialization ========== */
void Subband_InitAnalysisState(Subband_AnaState *st) {
    int k, t;
    for (k = 0; k < SUBBAND_D; k++) {
        for (t = 0; t < SUBBAND_L; t++) st->delay[k][t] = 0.0f;
        st->head[k] = 0;
    }
}
void Subband_InitSynthesisState(Subband_SynState *st) {
    int k, t;
    for (k = 0; k < SUBBAND_D; k++) {
        for (t = 0; t < SUBBAND_L; t++) st->delay[k][t] = 0.0f;
        st->head[k] = 0;
    }
}

/* ========== 8-point DFT (LUT-based) ========== */
static void dft8(const float re_in[], const float im_in[], float re_out[], float im_out[]) {
    int k, n;
    for (k = 0; k < SUBBAND_D; k++) {
        float sr = 0, si = 0;
        for (n = 0; n < SUBBAND_D; n++) {
            float c = W8_COS[k][n], s = W8_SIN[k][n];
            sr += re_in[n] * c + im_in[n] * s;    /* exp(-j*2*pi*k*n/D) */
            si += im_in[n] * c - re_in[n] * s;
        }
        re_out[k] = sr; im_out[k] = si;
    }
}

/* ========== 8-point IDFT (LUT-based) ========== */
static void idft8(const float re_in[], const float im_in[], float re_out[], float im_out[]) {
    int k, n;
    for (n = 0; n < SUBBAND_D; n++) {
        float sr = 0, si = 0;
        for (k = 0; k < SUBBAND_D; k++) {
            float c = W8_COS[k][n], s = W8_SIN[k][n];
            sr += re_in[k] * c - im_in[k] * s;    /* exp(+j*2*pi*k*n/D) */
            si += re_in[k] * s + im_in[k] * c;
        }
        re_out[n] = sr / SUBBAND_D; im_out[n] = si / SUBBAND_D;
    }
}

/* ========== Analysis Block: D inputs -> D complex sub-band outputs ========== */
void Subband_AnalysisBlock(Subband_AnaState *st, const int16_t in[SUBBAND_D], float out[SUBBAND_D][2]) {
    int k, tap;
    float u_re[SUBBAND_D], u_im[SUBBAND_D];
    float dft_re[SUBBAND_D], dft_im[SUBBAND_D];

    for (k = 0; k < SUBBAND_D; k++) {
        /* write input sample, advance write pointer (ring buffer) */
        st->delay[k][st->head[k]] = (float)in[k];
        st->head[k] = (st->head[k] + 1) & (SUBBAND_L - 1);

        /* L-tap polyphase convolution on branch k */
        float acc = 0;
        for (tap = 0; tap < SUBBAND_L; tap++) {
            int idx = (st->head[k] + SUBBAND_L - 1 - tap) & (SUBBAND_L - 1);
            acc += ana_poly[k][tap] * st->delay[k][idx];
        }
        u_re[k] = acc;
        u_im[k] = 0;
    }

    dft8(u_re, u_im, dft_re, dft_im);

    for (k = 0; k < SUBBAND_D; k++) {
        out[k][0] = dft_re[k]; out[k][1] = dft_im[k];
    }
}

/* ========== Synthesis Block: D complex sub-bands -> D outputs ========== */
void Subband_SynthesisBlock(Subband_SynState *st, const float in[SUBBAND_D][2], int16_t out[SUBBAND_D]) {
    int k, tap;
    float w_re[SUBBAND_D], w_im[SUBBAND_D];
    float idft_re[SUBBAND_D], idft_im[SUBBAND_D];

    for (k = 0; k < SUBBAND_D; k++) {
        w_re[k] = in[k][0]; w_im[k] = in[k][1];
    }

    idft8(w_re, w_im, idft_re, idft_im);

    for (k = 0; k < SUBBAND_D; k++) {
        /* write IDFT output, advance write pointer */
        st->delay[k][st->head[k]] = idft_re[k];
        st->head[k] = (st->head[k] + 1) & (SUBBAND_L - 1);

        /* L-tap polyphase convolution on branch k */
        float acc = 0;
        for (tap = 0; tap < SUBBAND_L; tap++) {
            int idx = (st->head[k] + SUBBAND_L - 1 - tap) & (SUBBAND_L - 1);
            acc += syn_poly[k][tap] * st->delay[k][idx];
        }

        /* saturate to int16 range */
        if (acc >  32767.0f) acc =  32767.0f;
        if (acc < -32768.0f) acc = -32768.0f;
        out[k] = (int16_t)acc;
    }
}

/* ========== Frame-level Pipeline: analyze -> passthrough -> synthesize ========== */
void Subband_Pipeline(Subband_AnaState *ana, Subband_SynState *syn,
                      const int16_t input[], int16_t output[], int len) {
    int blk_cnt = len / SUBBAND_D;
    int k, blk;
    float sbuf[SUBBAND_D][2];

    if (blk_cnt > SUBBAND_MAX_BLKS) blk_cnt = SUBBAND_MAX_BLKS;

    /* forward path: analysis */
    for (blk = 0; blk < blk_cnt; blk++) {
        Subband_AnalysisBlock(ana, &input[blk * SUBBAND_D], sbuf);
        for (k = 0; k < SUBBAND_D; k++) {
            pipe_buf[k][blk][0] = sbuf[k][0];
            pipe_buf[k][blk][1] = sbuf[k][1];
        }
    }
    /* backward path: synthesis (sub-band processing can be inserted here) */
    for (blk = 0; blk < blk_cnt; blk++) {
        for (k = 0; k < SUBBAND_D; k++) {
            sbuf[k][0] = pipe_buf[k][blk][0];
            sbuf[k][1] = pipe_buf[k][blk][1];
        }
        Subband_SynthesisBlock(syn, sbuf, &output[blk * SUBBAND_D]);
    }
}

/* ========== Main Example: ADC -> filter bank -> DAC, with key/LED control ========== */
void Subband_Example(void) {
    unsigned char ad_done = 0;
    static Subband_AnaState ana;
    static Subband_SynState syn;
    static int16_t adc_buf[2048], dac_buf[2048];

    Sys_Init();
    Led_Init();
    Key_Init();
    Led_Control(LED1_CORE, LED_ON);
    Led_Control(LED2_CORE, LED_ON);

    Subband_InitFilterCoeffs();
    Subband_InitAnalysisState(&ana);
    Subband_InitSynthesisState(&syn);

    Adc_Init(ADC_50KHZ, 2048);
    Dac_Init(DAC_50KHZ, 2048, DAC_CHANNEL_12);
    Adc_Start();
    Dac_Start();

    while (1) {
        if (FLAG_AD == 1) {
            FLAG_AD = 0;
            ad_done = 1;

            if (AD_Ping_Pong == AD_BUFFER_PONG)
                memcpy(adc_buf, AD_CH1_Buf0, sizeof(int16_t) * 2048);
            else
                memcpy(adc_buf, AD_CH1_Buf1, sizeof(int16_t) * 2048);
        }

        if (FLAG_DA == 1 && ad_done == 1) {
            FLAG_DA = 0;
            ad_done = 0;

            Subband_Pipeline(&ana, &syn, adc_buf, dac_buf, 2048);

            if (DA_Ping_Pong == DA_BUFFER_PONG) {
                memcpy(DA_CH1_Buf0, dac_buf, sizeof(int16_t) * 2048);
                memcpy(DA_CH2_Buf0, dac_buf, sizeof(int16_t) * 2048);
            } else {
                memcpy(DA_CH1_Buf1, dac_buf, sizeof(int16_t) * 2048);
                memcpy(DA_CH2_Buf1, dac_buf, sizeof(int16_t) * 2048);
            }

            Led_Control(LED1_BLUE, LED_TOGGLE);
        }

        if (FLAG_KEY1 == 1) { FLAG_KEY1 = 0; Adc_Start(); }
        if (FLAG_KEY2 == 1) { FLAG_KEY2 = 0; Adc_Stop();  }
        if (FLAG_KEY3 == 1) { FLAG_KEY3 = 0; Dac_Start(); }
        if (FLAG_KEY4 == 1) {
            FLAG_KEY4 = 0;
            Dac_Stop(); Adc_Stop();
            Subband_InitAnalysisState(&ana);
            Subband_InitSynthesisState(&syn);
            Adc_Start(); Dac_Start();
        }
        if (FLAG_KEY5 == 1) {
            FLAG_KEY5 = 0;
            Led_Control(LED1_RED, LED_TOGGLE);
            Led_Control(LED2_RED, LED_TOGGLE);
        }
    }
}