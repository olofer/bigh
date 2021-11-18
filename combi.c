/*

  combi.c

  Test program for big.h 

  Evaluate combinatorial numbers that tend to not fit in atomic integer types.
  Also evaluates basic integer operations: mult, div, mod, gcd, pow.
  Input arguments and print results in arbitrary digit bases (2..36).

  BUILD: gcc -Wall -O2 -o combi.exe combi.c

  USAGE: combi key1=val1 key2=val2 ...

  KEYS:
    inbase   : digit base for input numeric text (default=10)
    outbase  : digit base for output numeric text (default=10)
    words    : how many 32 bit words to use (default=32, 1024 bits)

    op       : bell (n)
               stirling (n, k)
               nchoosek (n, k)
               factorial (n)
               mult (n * k)
               gcd (n, k)
               mod (n % k)
               div (n / k)
               pow (n ^ k)

    n        : numeric text argument 1
    k        : numeric text argument 2 (if needed)

  EXAMPLES:
    combi n=100 k=50 op=nchoosek
    combi n=54 op=bell
    combi op=stirling n=10 k=4

 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include "big.h"

bool split_key_eq_value(const char* s, 
                        char* key, 
                        char* val)
{
  const char* q = s;
  while (true) {
    if (*q == '=' || *q == '\0') break;
    *key++ = *q++;
  }
  *key = '\0';
  if (*q != '=')
    return false;
  q++;
  while (*q != '\0') {
    *val++ = *q++;
  }
  *val = '\0';
  return true;
}

int main(int argc, char** argv)
{
  if (argc == 1) {
    printf("USAGE: %s key1=arg1 ...\n", argv[0]);
    return 0;
  }

  const int default_words = 32;
  const int default_inbase = 10;
  const int default_outbase = 10;

  const int strbufsz = 1024;
  char strbuf1[strbufsz];
  char strbuf2[strbufsz];
  char opstr[strbufsz];

  opstr[0] = '\0';

  int words = default_words;
  int inbase = default_inbase;
  int outbase = default_outbase;

  // Locate words and bases (key, value) pairs first since 
  // they are needed before reading actual numbers n, k
  for (int a = 1; a < argc; a++) {
    const bool oksplit = split_key_eq_value(argv[a], strbuf1, strbuf2);
    if (!oksplit) continue;
    if (strcmp(strbuf1, "words") == 0) {
      words = (int) strtol(strbuf2, NULL, 0);
    }
    if (strcmp(strbuf1, "inbase") == 0) {
      inbase = (int) strtol(strbuf2, NULL, 0);
    }
    if (strcmp(strbuf1, "outbase") == 0) {
      outbase = (int) strtol(strbuf2, NULL, 0);
    }
    if (strcmp(strbuf1, "op") == 0) {
      strcpy(opstr, strbuf2);
    }
  }

  if (words <= 0 || inbase < 2 || inbase > 36 || outbase < 2 || outbase > 36) {
    printf("invalid numeric inputs\n");
    return 1;
  }

  if (strlen(opstr) == 0) {
    printf("no op provided\n");
    return 1;
  }

  uint32_t N[words];
  uint32_t K[words];

  // If ncode (kcode) remains nonzero, then N (K) was not provided successfully.
  // If n, k are specified several times, the first successful will be kept
  uint32_t ncode = 0xffff;
  uint32_t kcode = 0xffff;
  for (int a = 1; a < argc; a++) {
    const bool oksplit = split_key_eq_value(argv[a], strbuf1, strbuf2);
    if (!oksplit) continue;
    if (strcmp(strbuf1, "n") == 0) {
      if (ncode != 0)
        ncode = big_init_string(N, strbuf2, inbase, words);
    }
    if (strcmp(strbuf1, "k") == 0) {
      if (kcode != 0)
        kcode = big_init_string(K, strbuf2, inbase, words);
    }
  }

  uint32_t rcode = 0xffff;
  uint32_t R[words];
  uint32_t Q[words];

  if (strcmp(opstr, "factorial") == 0 && ncode == 0) {
    if (big_is_equal_word(N, N[0], words))
      rcode = big_factorial(R, N[0], words);

  } else if (strcmp(opstr, "bell") == 0 && ncode == 0) {
    if (big_is_equal_word(N, N[0], words))
      rcode = big_bell(R, N[0], words);

  } else if (strcmp(opstr, "stirling") == 0 && ncode == 0 && kcode == 0) {
    if (big_is_equal_word(N, N[0], words) && big_is_equal_word(K, K[0], words))
      rcode = big_stirling(R, N[0], K[0], words);

  } else if (strcmp(opstr, "nchoosek") == 0 && ncode == 0 && kcode == 0) {
    if (big_is_equal_word(N, N[0], words) && big_is_equal_word(K, K[0], words))
      rcode = big_nchoosek(R, N[0], K[0], words);

  } else if (strcmp(opstr, "mult") == 0 && ncode == 0 && kcode == 0) {
    rcode = big_mult_big(R, N, K, false, words);

  } else if (strcmp(opstr, "mod") == 0 && ncode == 0 && kcode == 0) {
    rcode = big_div_big(Q, R, N, K, words);

  } else if (strcmp(opstr, "div") == 0 && ncode == 0 && kcode == 0) {
    rcode = big_div_big(R, Q, N, K, words);

  } else if (strcmp(opstr, "gcd") == 0 && ncode == 0 && kcode == 0) {
    rcode = big_gcd(R, N, K, words);

  } else if (strcmp(opstr, "pow") == 0 && ncode == 0 && kcode == 0) {
    if (big_is_equal_word(N, N[0], words) && big_is_equal_word(K, K[0], words))
      rcode = big_pow_word(R, N[0], K[0], words);

  } else {
    printf("unrecognized op, or missing argument n, or k\n");
    return 1;
  }

  if (rcode != 0) {
    printf("results error (or input error)\n");
    return 1;
  }

  rcode = big_to_string(strbuf1, strbufsz, outbase, R, words);

  if (rcode != 0) {
    printf("results text rendition error\n");
    return 1;
  }

  printf("%s\n", strbuf1);
  return 0;
}
