/**
 * user_subband.h - 8-band Polyphase DFT Filter Bank
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

/* ---------- constants ---------- */
#define SUBBAND_D          8                   /* number of sub-bands            */
#define SUBBAND_ORDER      64                  /* prototype filter order         */
#define SUBBAND_L          (SUBBAND_ORDER / SUBBAND_D)  /* polyphase branch length */
#define SUBBAND_MAX_BLKS   256                 /* max blocks per pipeline frame  */

/* ---------- per-branch delay line state ---------- */
typedef struct {
    float delay[SUBBAND_D][SUBBAND_L];         /* circular delay buffer          */
    int   head[SUBBAND_D];                     /* write pointer per branch       */
} Subband_AnaState;

typedef struct {
    float delay[SUBBAND_D][SUBBAND_L];
    int   head[SUBBAND_D];
} Subband_SynState;

/* ---------- public API ---------- */
void Subband_InitFilterCoeffs(void);
void Subband_InitAnalysisState(Subband_AnaState *st);
void Subband_InitSynthesisState(Subband_SynState *st);
void Subband_Pipeline(Subband_AnaState *ana, Subband_SynState *syn,
                      const int16_t input[], int16_t output[], int len);
void Subband_Example(void);

#endif /* _USER_SUBBAND_H_ */