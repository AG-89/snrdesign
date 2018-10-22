/*
 * File: _coder_EH_iterate_api.h
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 22-Oct-2018 05:00:10
 */

#ifndef _CODER_EH_ITERATE_API_H
#define _CODER_EH_ITERATE_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_EH_iterate_api.h"

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
extern void EH_iterate(emxArray_real_T *E, emxArray_real_T *H, emxArray_real_T
  *buffer, real_T maxlag, real_T newrange_size, emxArray_real_T *Eout,
  emxArray_real_T *Hout);
extern void EH_iterate_api(const mxArray *prhs[5], const mxArray *plhs[2]);
extern void EH_iterate_atexit(void);
extern void EH_iterate_initialize(void);
extern void EH_iterate_terminate(void);
extern void EH_iterate_xil_terminate(void);

#endif

/*
 * File trailer for _coder_EH_iterate_api.h
 *
 * [EOF]
 */
