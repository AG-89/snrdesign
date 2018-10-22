/*
 * File: findpeaks_fast_emxutil.h
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 22-Oct-2018 04:23:25
 */

#ifndef FINDPEAKS_FAST_EMXUTIL_H
#define FINDPEAKS_FAST_EMXUTIL_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rtwtypes.h"
#include "findpeaks_fast_types.h"

/* Function Declarations */
extern void emxEnsureCapacity(emxArray__common *emxArray, int oldNumel, unsigned
  int elementSize);
extern void emxFree_real_T(emxArray_real_T **pEmxArray);
extern void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions);

#endif

/*
 * File trailer for findpeaks_fast_emxutil.h
 *
 * [EOF]
 */
