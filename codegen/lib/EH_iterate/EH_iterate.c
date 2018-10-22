/*
 * File: EH_iterate.c
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 22-Oct-2018 05:00:10
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "EH_iterate.h"
#include "EH_iterate_emxutil.h"

/* Function Definitions */

/*
 * perform iterative calculation for E & H
 * Arguments    : emxArray_real_T *E
 *                emxArray_real_T *H
 *                const emxArray_real_T *buffer
 *                double maxlag
 *                double newrange_size
 *                emxArray_real_T *Eout
 *                emxArray_real_T *Hout
 * Return Type  : void
 */
void EH_iterate(emxArray_real_T *E, emxArray_real_T *H, const emxArray_real_T
                *buffer, double maxlag, double newrange_size, emxArray_real_T
                *Eout, emxArray_real_T *Hout)
{
  int varargin_2;
  int buffer_len;
  int sN;
  int L;
  varargin_2 = buffer->size[1];
  buffer_len = buffer->size[1];
  for (sN = 0; sN < (int)((newrange_size - 1.0) + 1.0); sN++) {
    /* sample number in range */
    for (L = 0; L < (int)maxlag; L++) {
      /* lag */
      /* .* ones(1,maxlag); */
      /* EHC_newterm4 = EHC_newterm2; */
      E->data[L] = (E->data[L] + buffer->data[(int)((double)varargin_2 - (double)
        sN) - 1] * buffer->data[(int)((double)varargin_2 - (double)sN) - 1]) -
        buffer->data[(int)(((double)varargin_2 - 2.0 * (1.0 + (double)L)) -
                           (double)sN) - 1] * buffer->data[(int)(((double)
        varargin_2 - 2.0 * (1.0 + (double)L)) - (double)sN) - 1];
      H->data[L] = (H->data[L] + buffer->data[(int)((double)buffer_len - (double)
        sN) - 1] * buffer->data[(int)(((double)buffer_len - (1.0 + (double)L)) -
        (double)sN) - 1]) - buffer->data[(int)(((double)buffer_len - (1.0 +
        (double)L)) - (double)sN) - 1] * buffer->data[(int)(((double)buffer_len
        - 2.0 * (1.0 + (double)L)) - (double)sN) - 1];
    }
  }

  varargin_2 = Eout->size[0] * Eout->size[1];
  Eout->size[0] = 1;
  Eout->size[1] = E->size[1];
  emxEnsureCapacity((emxArray__common *)Eout, varargin_2, sizeof(double));
  buffer_len = E->size[0] * E->size[1];
  for (varargin_2 = 0; varargin_2 < buffer_len; varargin_2++) {
    Eout->data[varargin_2] = E->data[varargin_2];
  }

  varargin_2 = Hout->size[0] * Hout->size[1];
  Hout->size[0] = 1;
  Hout->size[1] = H->size[1];
  emxEnsureCapacity((emxArray__common *)Hout, varargin_2, sizeof(double));
  buffer_len = H->size[0] * H->size[1];
  for (varargin_2 = 0; varargin_2 < buffer_len; varargin_2++) {
    Hout->data[varargin_2] = H->data[varargin_2];
  }
}

/*
 * File trailer for EH_iterate.c
 *
 * [EOF]
 */
