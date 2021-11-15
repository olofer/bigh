/*

  sha2consts.c

  Test program for big.h 

  Fractional bits of square roots or cube roots of first k primes
  These are used internally in the SHA-2 cryptographic hash function.

  BUILD: gcc -Wall -O2 -o sha2consts.exe sha2consts.c

  USAGE: sha2consts [sqrt|cbrt] k base fbits

  EXAMPLES:
    sha2consts sqrt 8 16 32         (SHA-256 constants h)
    sha2consts cbrt 64 16 32        (SHA-256/224 constants k)
    sha2consts sqrt 16 16 64        (SHA-224 constants h: last half of 8 last numbers)
    sha2consts sqrt 8 16 64         (SHA-512 constants h)
    sha2consts cbrt 80 16 64        (SHA-512 constants k)

 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include "big.h"

void generate_first_k_primes(int* P, int k) {
  if (k <= 0) return;
  P[0] = 2;
  if (k == 1) return;
  P[1] = 3;
  int n = 2;
  for (int q = 3; n != k; q += 2) {
    int i = 1;
    for ( ; i < n; i++) {
      int a = q / P[i];
      if (a * P[i] == q) break;
    }
    if (i == n) {
      P[n++] = q;
    }
  }
}

int main(int argc, char** argv)
{
  if (argc != 5) {
    printf("USAGE: %s [sqrt|cbrt] k base fbits\n", argv[0]);
    return 0;
  }

  const int numwords = 32;
  const bool cube_root = (strcmp(argv[1], "cbrt") == 0);
  const int num_primes = (int) strtod(argv[2], NULL);
  const int base = (int) strtod(argv[3], NULL);
  const int fbits = (int) strtod(argv[4], NULL);

  if ((base < 2 || base > 36) ||
      (fbits < 1) || (num_primes < 1))
  {
    printf("unsupported argument(s)\n");
    return 1;
  }

  int P[num_primes];
  generate_first_k_primes(P, num_primes);

  for (int i = 0; i < num_primes; i++) {
    uint32_t p[numwords];
    uint32_t q[numwords];
    uint32_t r[numwords];
    big_init_word(p, (uint32_t) P[i], numwords);

    int ipart_bits = 0;
    int shift_bits = 0;

    if (cube_root) {
      big_icbrt_bisect(q, p, numwords);
      ipart_bits = big_bits_width(numwords) - big_bits_clz(q, numwords);
      shift_bits = 3 * fbits;
      big_shift_left(p, shift_bits, numwords);
      big_icbrt_bisect(r, p, numwords);
    } else {
      big_isqrt(q, p, numwords);
      ipart_bits = big_bits_width(numwords) - big_bits_clz(q, numwords);
      shift_bits = 2 * fbits;
      big_shift_left(p, shift_bits, numwords);
      big_isqrt(r, p, numwords);
    }

    if (ipart_bits + shift_bits > big_bits_width(numwords)) {
      printf("warning: register size = %i too small\n", big_bits_width(numwords));
    }

    big_shift_left(q, fbits, numwords);
    big_sub_big(r, q, numwords); // remove whole part (so that only fractional bits remain)

    const int strbufsz = 1024;
    char strbuf[strbufsz];
     
    big_to_string(strbuf, strbufsz, base, r, numwords);

    // If the string has initial zeros then prettify the printout with padding
    // but only bother for base=2 and base=16 when fbits is a multiple of 4
    int len = strlen(strbuf);
    char strpad[64];
    if (base == 2 && len < fbits) {
      int j = 0;
      while (j < fbits - len) strpad[j++] = '0';
      strpad[j] = '\0';
      printf("[%04i]: %06i: %s%s\n", i, P[i], strpad, strbuf);
    } else if (base == 16 && ((fbits >> 2) << 2 == fbits) && len < fbits >> 2) {
      int j = 0;
      while (j < (fbits >> 2) - len) strpad[j++] = '0';
      strpad[j] = '\0';
      printf("[%04i]: %06i: %s%s\n", i, P[i], strpad, strbuf);
    } else {
      printf("[%04i]: %06i: %s\n", i, P[i], strbuf);
    }
  }

  return 0;
}
