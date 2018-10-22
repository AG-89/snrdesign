/*
 * File: findpeaks_fast.c
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 22-Oct-2018 04:23:25
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "findpeaks_fast.h"
#include "findpeaks_fast_emxutil.h"

/* Function Definitions */

/*
 * custom findpeaks function much faster than stock
 * pks = y values, locs = x values (indices)
 * force inline
 * Arguments    : const emxArray_real_T *y
 *                emxArray_real_T *pks
 *                emxArray_real_T *locs
 * Return Type  : void
 */
void findpeaks_fast(const emxArray_real_T *y, emxArray_real_T *pks,
                    emxArray_real_T *locs)
{
  int i;
  int loop_ub;
  unsigned int num;
  i = pks->size[0] * pks->size[1];
  pks->size[0] = 1;
  pks->size[1] = y->size[1];
  emxEnsureCapacity((emxArray__common *)pks, i, sizeof(double));
  loop_ub = y->size[1];
  for (i = 0; i < loop_ub; i++) {
    pks->data[i] = 0.0;
  }

  /* allocate (faster than growing) */
  i = locs->size[0] * locs->size[1];
  locs->size[0] = 1;
  locs->size[1] = y->size[1];
  emxEnsureCapacity((emxArray__common *)locs, i, sizeof(double));
  loop_ub = y->size[1];
  for (i = 0; i < loop_ub; i++) {
    locs->data[i] = 0.0;
  }

  if (!(y->size[1] < 3)) {
    num = 0U;

    /* I couldn't find a way to make this faster */
    /* matches stock output well, may have extra peaks at mesas/end */
    for (i = 0; i <= y->size[1] - 3; i++) {
      /* exclude bounds like stock */
      if ((y->data[i + 1] - y->data[i] > 0.0) && (y->data[i + 2] - y->data[i + 1]
           <= 0.0)) {
        num++;
        pks->data[(int)num - 1] = y->data[i + 1];
        locs->data[(int)num - 1] = 2.0 + (double)i;
      }
    }

    i = pks->size[0] * pks->size[1];
    if (1 > (int)num) {
      pks->size[1] = 0;
    } else {
      pks->size[1] = (int)num;
    }

    emxEnsureCapacity((emxArray__common *)pks, i, sizeof(double));

    /* output right sized array */
    i = locs->size[0] * locs->size[1];
    if (1 > (int)num) {
      locs->size[1] = 0;
    } else {
      locs->size[1] = (int)num;
    }

    emxEnsureCapacity((emxArray__common *)locs, i, sizeof(double));
  } else {
    /* require at least 3 points */
  }
}

/*
 * File trailer for findpeaks_fast.c
 *
 * [EOF]
 */
