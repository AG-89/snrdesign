/*
 * File: main.c
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 22-Oct-2018 05:00:10
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/
/* Include Files */
#include "rt_nonfinite.h"
#include "EH_iterate.h"
#include "main.h"
#include "EH_iterate_terminate.h"
#include "EH_iterate_emxAPI.h"
#include "EH_iterate_initialize.h"

/* Function Declarations */
static emxArray_real_T *argInit_1xUnbounded_real_T(void);
static double argInit_real_T(void);
static void main_EH_iterate(void);

/* Function Definitions */

/*
 * Arguments    : void
 * Return Type  : emxArray_real_T *
 */
static emxArray_real_T *argInit_1xUnbounded_real_T(void)
{
  emxArray_real_T *result;
  static int iv0[2] = { 1, 2 };

  int idx1;

  /* Set the size of the array.
     Change this size to the value that the application requires. */
  result = emxCreateND_real_T(2, iv0);

  /* Loop over the array to initialize each element. */
  for (idx1 = 0; idx1 < result->size[1U]; idx1++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result->data[result->size[0] * idx1] = argInit_real_T();
  }

  return result;
}

/*
 * Arguments    : void
 * Return Type  : double
 */
static double argInit_real_T(void)
{
  return 0.0;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
static void main_EH_iterate(void)
{
  emxArray_real_T *Eout;
  emxArray_real_T *Hout;
  emxArray_real_T *E;
  emxArray_real_T *H;
  emxArray_real_T *buffer;
  emxInitArray_real_T(&Eout, 2);
  emxInitArray_real_T(&Hout, 2);

  /* Initialize function 'EH_iterate' input arguments. */
  /* Initialize function input argument 'E'. */
  E = argInit_1xUnbounded_real_T();

  /* Initialize function input argument 'H'. */
  H = argInit_1xUnbounded_real_T();

  /* Initialize function input argument 'buffer'. */
  buffer = argInit_1xUnbounded_real_T();

  /* Call the entry-point 'EH_iterate'. */
  EH_iterate(E, H, buffer, argInit_real_T(), argInit_real_T(), Eout, Hout);
  emxDestroyArray_real_T(Hout);
  emxDestroyArray_real_T(Eout);
  emxDestroyArray_real_T(buffer);
  emxDestroyArray_real_T(H);
  emxDestroyArray_real_T(E);
}

/*
 * Arguments    : int argc
 *                const char * const argv[]
 * Return Type  : int
 */
int main(int argc, const char * const argv[])
{
  (void)argc;
  (void)argv;

  /* Initialize the application.
     You do not need to do this more than one time. */
  EH_iterate_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_EH_iterate();

  /* Terminate the application.
     You do not need to do this more than one time. */
  EH_iterate_terminate();
  return 0;
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
