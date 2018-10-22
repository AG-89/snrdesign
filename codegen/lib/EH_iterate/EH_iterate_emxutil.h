/*
 * File: EH_iterate_emxutil.h
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 22-Oct-2018 05:00:10
 */

#ifndef EH_ITERATE_EMXUTIL_H
#define EH_ITERATE_EMXUTIL_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rtwtypes.h"
#include "EH_iterate_types.h"

/* Function Declarations */
extern void emxEnsureCapacity(emxArray__common *emxArray, int oldNumel, unsigned
  int elementSize);
extern void emxFree_real_T(emxArray_real_T **pEmxArray);
extern void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions);

#endif

/*
 * File trailer for EH_iterate_emxutil.h
 *
 * [EOF]
 */
