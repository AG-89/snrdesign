/*
 * File: _coder_findpeaks_fast_info.c
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 22-Oct-2018 04:23:25
 */

/* Include Files */
#include "_coder_findpeaks_fast_info.h"

/* Function Definitions */

/*
 * Arguments    : void
 * Return Type  : const mxArray *
 */
const mxArray *emlrtMexFcnResolvedFunctionsInfo(void)
{
  const mxArray *nameCaptureInfo;
  const char * data[5] = {
    "789ce554cb4ac340149d6a7d2d14170a2efc028b9de2037c6cec130bad148d228868d2dcdae03c4a666a6b37fa217e804b7fc68d7fe24ad324d3a6852145d08d"
    "17c2993327997b72b8094a94ab0984d0120a6a792ec045c5439c42a335ae2734a86a0625479e53fa4b8875ce2474654088c3e0a44d2d703dc24c0a83636c4e1d",
    "6632693cb600b920387900db571a0e01c3a150e11172ec78849622d280f4a5fe3adf84fafd599b22b72986764994a0483e6f9af74f4e984f45934fffbeafa7a1"
    "7e55bcce1fe07301aec0794e2d2f118b5b4c6e6db6bd3d5ce0f5360526056e7670315fdcd9decf6470356b54b239dc7098dd02f35edc344c21d334eabfabf137",
    "3ba1ff698dff8510d74676578e025c0f71e328cec7a439ae6a7ca81c955ee736b869c71b2e9799244d80ddc9261ae671fb431fe3a5f3a14af57bfd613f757e2e"
    "a69fd2afca954b7f84b2427057828daba624a68525e7c4e25d0c94f857ca4f08a75444388828989af89c269d1bdd7f6101cd87abe7c3cf8ff7ecdff50beabff4",
    "fbabefad5aead42edc9a51de2dec9d5a27bbbd96dd93111fb5983e713e9086fff6f9df8d3482e8",
    "" };

  nameCaptureInfo = NULL;
  emlrtNameCaptureMxArrayR2016a(data, 1832U, &nameCaptureInfo);
  return nameCaptureInfo;
}

/*
 * Arguments    : void
 * Return Type  : mxArray *
 */
mxArray *emlrtMexFcnProperties(void)
{
  mxArray *xResult;
  mxArray *xEntryPoints;
  const char * fldNames[4] = { "Name", "NumberOfInputs", "NumberOfOutputs",
    "ConstantInputs" };

  mxArray *xInputs;
  const char * b_fldNames[4] = { "Version", "ResolvedFunctions", "EntryPoints",
    "CoverageInfo" };

  xEntryPoints = emlrtCreateStructMatrix(1, 1, 4, fldNames);
  xInputs = emlrtCreateLogicalMatrix(1, 1);
  emlrtSetField(xEntryPoints, 0, "Name", mxCreateString("findpeaks_fast"));
  emlrtSetField(xEntryPoints, 0, "NumberOfInputs", mxCreateDoubleScalar(1.0));
  emlrtSetField(xEntryPoints, 0, "NumberOfOutputs", mxCreateDoubleScalar(2.0));
  emlrtSetField(xEntryPoints, 0, "ConstantInputs", xInputs);
  xResult = emlrtCreateStructMatrix(1, 1, 4, b_fldNames);
  emlrtSetField(xResult, 0, "Version", mxCreateString("9.2.0.538062 (R2017a)"));
  emlrtSetField(xResult, 0, "ResolvedFunctions", (mxArray *)
                emlrtMexFcnResolvedFunctionsInfo());
  emlrtSetField(xResult, 0, "EntryPoints", xEntryPoints);
  return xResult;
}

/*
 * File trailer for _coder_findpeaks_fast_info.c
 *
 * [EOF]
 */
