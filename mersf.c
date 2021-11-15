/*
 * mersf.c
 *
 * Demo / test of big.h codebase.
 * 
 * Factorize a Mersenne number using big integer 
 * arithmetic (32*w bits) and present factors
 * with base-b digits.
 * Or factorize a general integer provided as a 
 * base-b string input n.
 *
 * M(n) = 2^n - 1, specify integer n.
 * M(n) itself is returned as a base-b string x.
 * Then factorize M(n) and return all prime factors
 * as base-b digit strings in cell array y.
 * Optionally test z = product of y.
 *
 * USAGE:
 *   [x, y, z] = mersf(n, b, w);
 *
 * BUILD:
 *   MATLAB: mex mersf.c
 *   OCTAVE: mkoctfile --mex --strip mersf.c
 *
 * BACKGROUND:
 *   https://en.wikipedia.org/wiki/Mersenne_prime
 *
 */

#include "mex.h"

#include <stdint.h>
#include "big.h"

#define THEPRINTF mexPrintf
#define THEERRMSG mexErrMsgTxt
#define THEWRNMSG mexWarnMsgTxt

/* ------------------------------------------------------------ */
/* Some MEX I/O helpers                                         */
/* ------------------------------------------------------------ */

#define MIN_ARGUMENTS 1
#define MAX_ARGUMENTS 3
#define ARG_N prhs[0]
#define ARG_B prhs[1]
#define ARG_W prhs[2]

#define MAX_OUTARGS 3
#define OUT_X plhs[0]
#define OUT_Y plhs[1]
#define OUT_Z plhs[2]

/* ------------------------------------------------------------ */
/* MEX main                                                     */
/* ------------------------------------------------------------ */

void mexFunction(
  int nlhs, mxArray *plhs[],
  int nrhs, const mxArray *prhs[])
{
  const int default_base = 10;
  const int default_numwords = 32;

  if (nrhs < MIN_ARGUMENTS || nrhs > MAX_ARGUMENTS) {
    THEPRINTF("USAGE: [x, y, z] = mersf(n [, b, w]);\n");
    THEPRINTF("PURPOSE: factorize x = M(n) = 2^n - 1, and show base-b digits.\n");
    THEPRINTF("For general integer factorization: set n to base-b string.\n");
    THEPRINTF("Default base b = %i, default register size w = %i (%i bits).\n", 
      default_base, default_numwords, big_bits_width(default_numwords));
    #ifdef HAVE_OCTAVE
    THEPRINTF("(compiled in Octave)\n");
    #else
    THEPRINTF("(compiled in Matlab)\n");
    #endif
    THEERRMSG("Incorrect number of input arguments provided");
  }

  if (nlhs > MAX_OUTARGS) {
    THEERRMSG("Too many output argument specified");
  }

  int base = default_base;
  if (nrhs >= 2) {
    int db = (int) *mxGetPr(ARG_B);
    if (db > 1) base = db;
    if (base != default_base) {
      THEPRINTF("setting general digit base = %i\n", base);  
    }
  }

  int numwords = default_numwords;

  if (nrhs == 3) {
    int nw = (int) *mxGetPr(ARG_W);
    if (nw > 0) numwords = nw;
    if (numwords != default_numwords) {
      THEPRINTF("setting register width = %i bits (%i words)\n", big_bits_width(numwords), numwords);
    }
  }

  uint32_t M[numwords];

  const char* strN = mxArrayToString(ARG_N);
  if (strN != NULL) {
    // String input n; try to interpret as base-b integer
    const uint32_t okN = big_init_string(M, strN, base, numwords);
    if (okN != 0) {
      THEERRMSG("Failed to get numerical representation from string n.");
    }
  } else {
    // Input a number n; generate Mersenne number M(n)
    const double* pN = mxGetPr(ARG_N);
    const int nN = mxGetNumberOfElements(ARG_N);
    if (nN != 1 || pN == NULL) {
      THEERRMSG("scalar n required\n");
    }
    const int n = (int) *pN;
    if (n < 2) {
      THEERRMSG("n >= 2 required\n");
    }
    big_init_word(M, (uint32_t) 1, numwords);
    big_shift_left(M, n, numwords);
    big_sub_word(M, (uint32_t) 1, numwords);
  }

  const int strbufsz = 10000;
  char strbuf[strbufsz];

  uint32_t rval = big_to_string(strbuf, strbufsz, base, M, numwords);
  OUT_X = mxCreateString(strbuf);

  if (rval != 0) {
    THEERRMSG("Failed to create string representation of big integer.\n");
  }

  if (nlhs < 2) return;

  const int max_factors = 2000;
  int num_factors = 0;
  mxArray* factorStrings[max_factors];

  uint32_t Q[numwords];
  uint32_t R[numwords];
  uint32_t N[numwords];
  uint32_t sqrtM[numwords];
  big_isqrt_bisect(sqrtM, M, numwords);
  big_init_word(N, 0, numwords);
  while (big_is_less_than(N, sqrtM, true, numwords)) {
    if (big_is_equal_word(N, 2, numwords))
      big_add_word(N, 1, numwords);
    else
      big_add_word(N, 2, numwords);
    // the sequence of N will always be 2, 3, 5, 7, ...
    big_div_big(Q, R, M, N, numwords);
    if (big_is_zero(R, numwords)) {
      big_to_string(strbuf, strbufsz, base, N, numwords);
      factorStrings[num_factors++] = mxCreateString(strbuf);
      big_init_big(M, Q, numwords);
      big_isqrt_bisect(sqrtM, M, numwords);
      big_init_word(N, 0, numwords);
    }
  }
  // if the largest factor is 2, we end up with M=1 here; which should be ignored
  if (!big_is_equal_word(M, 1, numwords)) {
    big_to_string(strbuf, strbufsz, base, M, numwords);
    factorStrings[num_factors++] = mxCreateString(strbuf);
  }
  // Create cell array of strings assigned to output
  OUT_Y = mxCreateCellMatrix(num_factors, 1);
  for (int i = 0; i < num_factors; i++) {
    mxSetCell(OUT_Y, i, factorStrings[i]);
  }

  if (nlhs < 3) return;

  // if the 3rd output is requested; re-read all string factors and multiply them
  big_init_string(M, mxArrayToString(mxGetCell(OUT_Y, 0)), base, numwords);
  for (int i = 1; i < num_factors; i++) {
    big_init_string(N, mxArrayToString(mxGetCell(OUT_Y, i)), base, numwords);
    big_mult(M, N, numwords);
  }

  big_to_string(strbuf, strbufsz, base, M, numwords);
  OUT_Z = mxCreateString(strbuf);

  return;
}
