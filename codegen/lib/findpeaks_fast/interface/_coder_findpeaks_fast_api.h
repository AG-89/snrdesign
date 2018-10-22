/*
 * File: _coder_findpeaks_fast_api.h
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 22-Oct-2018 04:23:25
 */

#ifndef _CODER_FINDPEAKS_FAST_API_H
#define _CODER_FINDPEAKS_FAST_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_findpeaks_fast_api.h"

/* Type Definitions */
#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T

struct emxArray_real_T
{
  real_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray_real_T*/

#ifndef typedef_emxArray_real_T
#define typedef_emxArray_real_T

typedef struct emxArray_real_T emxArray_real_T;

#endif                                 /*typedef_emxArray_real_T*/

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void findpeaks_fast(emxArray_real_T *y, emxArray_real_T *pks,
  emxArray_real_T *locs);
extern void findpeaks_fast_api(const mxArray * const prhs[1], const mxArray
  *plhs[2]);
extern void findpeaks_fast_atexit(void);
extern void findpeaks_fast_initialize(void);
extern void findpeaks_fast_terminate(void);
extern void findpeaks_fast_xil_terminate(void);

#endif

/*
 * File trailer for _coder_findpeaks_fast_api.h
 *
 * [EOF]
 */
