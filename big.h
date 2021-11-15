#ifndef __BIG_H__
#define __BIG_H__

/* ------------------------------------------------------------ */
/* Basic header-only support for big integer arithmetic.        */
/* Maybe educational, but performance not guaranteed.           */
/*                                                              */
/* Requires <stdint.h>, <stdbool.h>                             */
/* ------------------------------------------------------------ */

void big_init_word(uint32_t* x, uint32_t y, int w) {
  x[0] = y;
  for (int i = 1; i < w; i++)
    x[i] = 0;
}

void big_init_shifted_word(uint32_t* x, uint32_t y, int s, int w) {
  for (int i = 0; i < w; i++)
    x[i] = 0;
  if (s >= 0 && s < w)
    x[s] = y;
}

void big_init_long(uint32_t * x, uint64_t y, int w) {
  x[0] = (uint32_t) (y & 0xffffffff);
  x[1] = (uint32_t) (y >> 32);
  for (int i = 2; i < w; i++)
    x[i] = 0;
}

uint32_t big_init_double_(uint32_t* x, double y, int w) {
  // TODO: direct conversion from bit patterns
  // https://en.wikipedia.org/wiki/Double-precision_floating-point_format#IEEE_754_double-precision_binary_floating-point_format:_binary64
  return 1;
}

// clumsy but ~works; more efficient (but less portable)
// to directly extract the mantissa and exponent bits from the double
uint32_t big_init_double(uint32_t* x, double y, int w) {
  big_init_word(x, 0, w);
  const double two_pow_32 = (double)(((uint64_t) 1) << 32);
  if (y < 0) return 0;
  double y1 = two_pow_32;
  int idx = 0;
  for ( ; idx < w; idx++) {
    if (y < y1) break;
    y1 *= two_pow_32;
  }
  if (idx == w) return idx;
  y1 /= two_pow_32;
  while (idx >= 0) {
    double q = y / y1;
    uint32_t d = (uint32_t) q;
    x[idx--] = d;
    y -= d * y1;
    y1 /= two_pow_32;
  }
  return (idx + 1);
}

void big_init_big(uint32_t* x, const uint32_t* y, int w) {
  for (int i = 0; i < w; i++)
    x[i] = y[i];
}

// sometimes useful, but lossy and overflows if pushed too far
// (more efficient to directly generate IEEE bit pattern)
double big_to_double(const uint32_t* x, int w) {
  const double two_pow_32 = (double)(((uint64_t) 1) << 32);
  double s = 1.0;
  double y = (double) x[0];
  for (int i = 1; i < w; i++) {
    s *= two_pow_32;
    y += s * (double) x[i];
  }
  return y;
}

int big_most_significant_word(const uint32_t* x, int w) {
  for (int i = w - 1; i >= 0; i--)
    if (x[i] != 0)
      return i;
  return -1; // if x is zero, -1 is returned
}

bool big_is_equal(const uint32_t* x, 
                  const uint32_t* y, 
                  int w) 
{
  for (int i = 0; i < w; i++)
    if (x[i] != y[i])
      return false;
  return true;
}

bool big_is_equal_word(const uint32_t* x, 
                       uint32_t y, 
                       int w) 
{
  uint32_t ybig[w];
  big_init_word(ybig, y, w);
  return big_is_equal(x, ybig, w);
}

bool big_is_greater_than(const uint32_t* x, 
                         const uint32_t* y, 
                         bool or_equal, 
                         int w) 
{
  for (int i = w - 1; i >= 0; i--)
    if (x[i] != y[i])
      return (x[i] > y[i]);
  return or_equal;
}

bool big_is_less_than(const uint32_t* x, 
                      const uint32_t* y, 
                      bool or_equal, 
                      int w) 
{
  for (int i = w - 1; i >= 0; i--)
    if (x[i] != y[i])
      return (x[i] < y[i]);
  return or_equal;
}

bool big_is_zero(const uint32_t* x, int w) {
  for (int i = 0; i < w; i++)
    if (x[i] != 0) return false;
  return true;
}

void big_bits_and(uint32_t* x, const uint32_t* y, int w) {
  for (int i = 0; i < w; i++)
    x[i] &= y[i];
}

void big_bits_or(uint32_t* x, const uint32_t* y, int w) {
  for (int i = 0; i < w; i++)
    x[i] |= y[i];
}

void big_bits_xor(uint32_t* x, const uint32_t* y, int w) {
  for (int i = 0; i < w; i++)
    x[i] ^= y[i];
}

void big_bits_not(uint32_t* x, int w) {
  for (int i = 0; i < w; i++)
    x[i] = ~x[i];
}

uint32_t big_bits_width(int w) {
  return (8 * w * sizeof(uint32_t));
}

#if defined(__GNUC__)

uint32_t big_bits_count(const uint32_t* x, int w) {
  uint32_t c = 0;
  for (int i = 0; i < w; i++) {
    //uint32_t v = x[i];
    //for ( ; v; c++) v &= (v - 1);
    c += __builtin_popcount(x[i]);
  }
  return c;
}

uint32_t big_bits_clz(const uint32_t* x, int w) {
  uint32_t c = 0;
  for (int i = w - 1; i >= 0; i--) {
    if (x[i] == 0) {
      c += 32;
      continue;
    }
    c += __builtin_clz(x[i]);
    break; 
  }
  return c;
}

uint32_t big_bits_ctz(const uint32_t* x, int w) {
  uint32_t c = 0;
  for (int i = 0; i < w; i++) {
    if (x[i] == 0) {
      c += 32;
      continue;
    }
    c += __builtin_ctz(x[i]);
    break;
  }
  return c;
}

#endif

// multiply by 2^32, plus add in
uint32_t big_shift_word_up(uint32_t* x, uint32_t in, int w) {
  const uint32_t out = x[w - 1];
  for (int i = w - 1; i > 0; i--)
    x[i] = x[i - 1];
  x[0] = in;
  return out;
}

// divide by 2^32, if in = 0
uint32_t big_shift_word_down(uint32_t* x, uint32_t in, int w) {
  const uint32_t out = x[0];
  for (int i = 0; i < w - 1; i++)
    x[i] = x[i + 1];
  x[w - 1] = in;
  return out;
}

// y += a, a is just a single word
uint32_t big_add_word(uint32_t* y, 
                      uint32_t a, 
                      int w) 
{
  uint32_t carry = a;
  for (int i = 0; i < w; i++) {
    uint64_t z = (uint64_t) y[i] + carry;
    y[i] = (uint32_t)(z & 0xffffffff); 
    carry = (uint32_t)(z >> 32);
  }
  return carry;
}

// x *= a
uint32_t big_mult_word(uint32_t* x, 
                       uint32_t a, 
                       int w) 
{
  uint32_t carry = 0;
  for (int i = 0; i < w; i++) {
    uint64_t z = (uint64_t) x[i];
    z *= a;
    z += carry;
    x[i] = (uint32_t)(z & 0xffffffff); 
    carry = (uint32_t)(z >> 32);
  }
  return carry;
}

// x += y + a
uint32_t big_add_big(uint32_t* x, 
                     const uint32_t* y, 
                     uint32_t a, 
                     int w) 
{
  uint32_t carry = a;
  for (int i = 0; i < w; i++) {
    uint64_t z = (uint64_t) carry;
    z += x[i];
    z += y[i];
    x[i] = (uint32_t)(z & 0xffffffff); 
    carry = (uint32_t)(z >> 32);
  }
  return carry;
}

// x -= y, using the two-complement's method
uint32_t big_sub_big(uint32_t* x, 
                     const uint32_t* y, 
                     int w) 
{
  uint32_t tmp[w];
  big_init_big(tmp, y, w);
  big_bits_not(tmp, w);
  return big_add_big(x, tmp, 1, w);
}

// x -= y, using the two-complement's method
uint32_t big_sub_word(uint32_t* x, 
                      uint32_t y, 
                      int w) 
{
  uint32_t tmp[w];
  big_init_word(tmp, y, w);
  big_bits_not(tmp, w);
  return big_add_big(x, tmp, 1, w);
}

// if accum: z += x * y, if not: z = x * y
// returns nonzero when the results overflow
uint32_t big_mult_big(uint32_t* z, 
                      const uint32_t* x, 
                      const uint32_t* y, 
                      bool accum, 
                      int w) 
{
  uint32_t overflow = 0;
  uint32_t tmp[w];
  if (!accum) big_init_word(z, 0, w);
  for (int i = 0; i < w; i++) {
    if (y[i] == 0) continue;
    big_init_big(tmp, x, w);
    overflow = big_mult_word(tmp, y[i], w);
    if (overflow) break;
    for (int j = 0; j < i; j++)
      big_shift_word_up(tmp, 0, w);
    overflow = big_add_big(z, tmp, 0, w);
    if (overflow) break;
  }
  return overflow;
}

// x *= y
uint32_t big_mult(uint32_t* x, const uint32_t* y, int w) {
  uint32_t tmp[w];
  big_init_big(tmp, x, w);
  return big_mult_big(x, tmp, y, false, w);
}

// x /= y, and return value = remainder
// do not do anything if y == 0
uint32_t big_div_word(uint32_t* x, uint32_t y, int w) {
  if (y == 0) return 0;
  int j = big_most_significant_word(x, w);
  if (j == -1) return y;
  uint32_t r = 0;
  for ( ; j >= 0; j--) {
    const uint32_t xj = x[j];
    const uint64_t aj = (((uint64_t) r) << 32) | ((uint64_t) xj);
    const uint32_t qj = aj / y;
    r = (uint32_t) (aj - qj * y);
    x[j] = qj;
  }
  return r;
}

// return the word digit d such that x = 2^(32*i)*d * y + r, 0 <= r 
// compute d using binary search (and return residual r)
uint32_t big_div_big_digit(uint32_t i, 
                           uint32_t* r,
                           const uint32_t* x, 
                           const uint32_t* y, 
                           int w)
{
  big_init_big(r, y, w);
  if (i >= w) return 0;
  uint32_t p[w];
  uint32_t d = 0;
  for (int s = 31; s >= 0; s--) {
    uint32_t dtest = d | (((uint32_t) 1) << s);
    big_init_big(p, y, w);
    big_mult_word(p, dtest, w);
    for (int j = 0; j < i; j++)
      big_shift_word_up(p, 0, w);
    if (big_is_less_than(p, x, true, w)) {
      d = dtest; 
      big_init_big(r, x, w);
      big_sub_big(r, p, w);  // r = x - p, p <= x
    }
  }
  return d;
}

// floor(x / y) = q
// y * q + r = x 
// find largest q such that r >= 0
uint32_t big_div_big(uint32_t* q, 
                     uint32_t* r, 
                     const uint32_t* x, 
                     const uint32_t* y, 
                     int w) 
{
  big_init_word(q, 0, w);
  big_init_word(r, 0, w);
  if (big_is_zero(y, w)) return 1;
  if (big_is_less_than(x, y, false, w)) {
    big_init_big(r, x, w);
    return 0;
  }
  // Now x >= y > 0 so division makes sense
  uint32_t* tmp = r;
  big_init_big(tmp, y, w);
  int msw = -1;
  while (big_is_greater_than(x, tmp, true, w)) {
    big_shift_word_up(tmp, 0, w);
    msw++;
  }
  uint32_t xpart[w];
  big_init_big(r, x, w);
  while (msw >= 0) {
    big_init_big(xpart, r, w);
    uint32_t dig = big_div_big_digit(msw, r, xpart, y, w);
    q[msw--] = dig;
  }
  return 0;
}

// meaningful range: 0 < s < 32
// for multiples of 32 use the full word shift
uint32_t big_shift_bit_up(uint32_t* x, int s, int w) {
  const uint32_t out = x[w - 1] >> (32 - s);
  for (int i = w - 1; i > 0; i--) {
    const uint32_t x1 = x[i];
    const uint32_t x0 = x[i - 1];
    x[i] = (x1 << s) | (x0 >> (32 - s));
  }
  x[0] <<= s;
  return out;
}

// meaningful range: 0 < s < 32
// for multiples of 32 use the full word shift
uint32_t big_shift_bit_down(uint32_t* x, int s, int w) {
  const uint32_t out = x[0] << (32 - s);
  for (int i = 0; i < w - 1; i++) {
    const uint32_t x0 = x[i];
    const uint32_t x1 = x[i + 1];
    x[i] = (x1 << (32 - s)) | (x0 >> s);
  }
  x[w - 1] >>= s;
  return out;
}

// General case, s >= 0
void big_shift_left(uint32_t* x, int s, int w) {
  while (s >= 32) {
    big_shift_word_up(x, 0, w);
    s -= 32;
  }
  if (s > 0)
    big_shift_bit_up(x, s, w);
}

// General case, s >= 0
void big_shift_right(uint32_t* x, int s, int w) {
  while (s >= 32) {
    big_shift_word_down(x, 0, w);
    s -= 32;
  }
  if (s > 0)
    big_shift_bit_down(x, s, w);
}

// Return the integer base-2 logarithm of x
// find the largest r such that 2^r <= x
uint32_t big_ilog2(const uint32_t* x, int w) {
  uint32_t y[w];
  big_init_big(y, x, w);
  uint32_t r = 0;
  big_shift_bit_down(y, 1, w);
  while (!big_is_zero(y, w)) {
    if (y[0] == 0) {
      r += 32;
      big_shift_word_down(y, 0, w);
      continue;
    }
    r++;
    big_shift_bit_down(y, 1, w);
  }
  return r;
}

// x <- n! (return nonzero if overflowing)
uint32_t big_factorial(uint32_t* x, 
                       uint32_t n, 
                       int w) 
{
  big_init_word(x, 1, w);
  if (n == 0) return 0;
  for (uint32_t i = 1; i <= n; i++) {
    const uint32_t carry = big_mult_word(x, i, w);
    if (carry != 0)
      return carry;
  }
  return 0;
}

// nonzero return value means failure (overflow)
uint32_t big_nchoosek(uint32_t* x, 
                      uint32_t n, 
                      uint32_t k, 
                      int w) 
{
  if (k > n) {
    big_init_word(x, 0, w);
    return 1;
  }
  if (k == 0 || k == n) {
    big_init_word(x, 1, w);
    return 0;
  }
  if (n - k < k) k = n - k;  // nchoosek(n, k) == nchoosek(n, n - k)
  big_init_word(x, 1, w);
  for (uint32_t i = 0; i < k; i++) {
    uint32_t c1 = big_mult_word(x, n - i, w);
    if (c1 != 0) return c1;
  }
  for (uint32_t j = 1; j <= k; j++) {
    uint32_t c1 = big_div_word(x, j, w);
    if (c1 != 0) return c1;
  }
  return 0;
}

// https://en.wikipedia.org/wiki/Euclidean_algorithm
void big_gcd(uint32_t* x,
             const uint32_t* a,
             const uint32_t* b,
             int w)
{
  if (big_is_equal(a, b, w)) {
    big_init_big(x, a, w);
    return;
  }
  uint32_t A[w];
  uint32_t B[w];
  if (big_is_greater_than(a, b, false, w)) {
    big_init_big(A, a, w);
    big_init_big(B, b, w);
  } else {
    big_init_big(A, b, w);
    big_init_big(B, a, w);
  }
  // Ensured A > B
  uint32_t Q[w];
  uint32_t R[w];
  while (!big_is_zero(B, w)) {
    big_div_big(Q, R, A, B, w);  // R = A mod B
    big_init_big(A, B, w);
    big_init_big(B, R, w);
  }
  big_init_big(x, A, w);
}

// least common multiple x = lcm(a, b) = (a * b) / gcd(a, b)
void big_lcm(uint32_t* x,
             const uint32_t* a,
             const uint32_t* b,
             int w)
{
  uint32_t gcdab[w];
  uint32_t ab[w];
  uint32_t rem[w];
  big_gcd(gcdab, a, b, w);
  big_mult_big(ab, a, b, false, w);
  big_div_big(x, rem, ab, gcdab, w);
}

void big_isqrt(uint32_t* y,
               const uint32_t* x,
               int w)
{
  uint32_t num[w];
  uint32_t res[w];
  uint32_t bit[w];
  uint32_t z[w];
  big_init_word(res, 0, w);
  big_init_word(bit, 1, w);
  big_init_big(num, x, w);
  big_shift_left(bit, big_bits_width(w) - 2, w);
  while (big_is_greater_than(bit, num, false, w))
    big_shift_right(bit, 2, w);
  while (!big_is_zero(bit, w)) {
    big_init_big(z, res, w);
    big_add_big(z, bit, 0, w);      // z = res + bit
    if (big_is_greater_than(num, z, true, w)) {
      big_sub_big(num, z, w);       // num -= z
      big_shift_right(res, 1, w);   // res >>= 1
      big_add_big(res, bit, 0, w);  // res += bit
    } else {
      big_shift_right(res, 1, w);
    }
    big_shift_right(bit, 2, w);
  }
  big_init_big(y, res, w);
}

// y <- integer square root of x, using bisection; not very efficient
uint32_t big_isqrt_bisect(uint32_t* y, 
                          const uint32_t* x, 
                          int w) 
{
  big_init_word(y, 0, w);
  if (big_is_zero(x, w)) return 0;
  uint32_t a[w];
  uint32_t b[w];
  uint32_t c[w];
  uint32_t s = big_ilog2(x, w) >> 1;
  big_init_word(a, 1, w);
  big_shift_left(a, s, w);
  big_init_big(b, a, w);
  big_shift_left(b, 1, w);
  big_init_big(c, a, w);
  big_add_big(c, b, 0, w);
  big_shift_right(c, 1, w);
  while (!big_is_equal(a, c, w)) {
    const uint32_t carry = big_mult_big(y, c, c, false, w);
    if (carry != 0) return carry;
    if (big_is_greater_than(y, x, false, w)) {
      big_init_big(b, c, w);
    } else {
      big_init_big(a, c, w);
    }
    big_init_big(c, a, w);
    big_add_big(c, b, 0, w);
    big_shift_right(c, 1, w);
  }
  big_init_big(y, a, w);
  return 0;
}

// y <- integer cube root of x, using bisection
uint32_t big_icbrt_bisect(uint32_t* y, 
                          const uint32_t* x, 
                          int w) 
{
  big_init_word(y, 0, w);
  if (big_is_zero(x, w)) return 0;
  uint32_t a[w];
  uint32_t b[w];
  uint32_t c[w];
  uint32_t s = big_ilog2(x, w) / 3;
  big_init_word(a, 1, w);
  big_shift_left(a, s, w);
  big_init_big(b, a, w);
  big_shift_left(b, 1, w);
  big_init_big(c, a, w);
  big_add_big(c, b, 0, w);
  big_shift_right(c, 1, w);
  while (!big_is_equal(a, c, w)) {
    uint32_t carry = big_mult_big(y, c, c, false, w);
    if (carry != 0) return carry;
    carry = big_mult(y, c, w);  // y *= c
    if (carry != 0) return carry;
    if (big_is_greater_than(y, x, false, w)) {
      big_init_big(b, c, w);
    } else {
      big_init_big(a, c, w);
    }
    big_init_big(c, a, w);
    big_add_big(c, b, 0, w);
    big_shift_right(c, 1, w);
  }
  big_init_big(y, a, w);
  return 0;
}

/*
   Allow digit base >= 2 and <= 36 (0,..,9,A,..,Z)
 */

uint32_t big_to_string(char* s, 
                       int maxs, 
                       uint32_t b, 
                       const uint32_t* x, 
                       int w) 
{
  if (s == NULL || maxs <= 2) return 1;
  s[0] = '\0';
  if (b < 2 || b > 36) return 1;
  uint32_t B[w];
  big_init_word(B, 1, w);
  int msd = 0;
  while (big_is_less_than(B, x, true, w)) {
    const uint32_t c2 = big_mult_word(B, b, w);  // B *= b
    if (c2 != 0) return c2;
    msd++;
  }
  if (msd == 0) {
    s[0] = (char) 48;
    s[1] = '\0';
    return 0;
  }
  uint32_t xpart[w];
  uint32_t r[w];
  big_init_big(r, x, w);
  uint32_t q[w];
  // Now b^msd is strictly larger than x
  if (msd >= maxs) return 1;  // string buffer is not big enough
  msd--;
  for (int i = 0; i <= msd; i++) {
    uint32_t c = big_div_word(B, b, w); // go digit lower
    if (c != 0) return 1;
    big_init_big(xpart, r, w);
    c = big_div_big(q, r, xpart, B, w);
    if (c != 0) return 1;
    uint8_t digit = (uint8_t) q[0];
    if (digit >= b) return 1;
    s[i] = (char) (digit < 10 ? (48 + digit) : (65 + digit - 10));
  }
  s[msd + 1] = '\0';
  return 0;
}

// Return nonzero when there is error; 
// b required to be >= 2 and <= 10 (could be extended)
uint32_t big_init_string(uint32_t* x, 
                         const char* s, 
                         uint32_t b, 
                         int w) 
{
  big_init_word(x, 0, w);
  if (s == NULL) return 1;
  if (b < 2 || b > 36) return 1;
  uint32_t B[w];
  big_init_word(B, 1, w);
  uint32_t y[w];
  int idx = 0;
  for ( ; s[idx] != '\0'; idx++);
  for (int i = idx - 1; i >= 0; i--) {
    uint8_t d = (uint8_t) s[i];
    if (d >= 48 && d <= 57) d -= 48;  // '0',...,'9': d = 0..9
      else if (d >= 65 && d <= 90) d = d - 65 + 10;  // 'A',...,'Z': d = 10..35
        else return 1; 
    if (d >= b) return 1;
    big_init_word(y, (uint32_t) d, w);
    const uint32_t c1 = big_mult_big(x, B, y, true, w);  // x += B * y
    if (c1 != 0) return c1;
    const uint32_t c2 = big_mult_word(B, b, w);          // B *= b
    if (c2 != 0) return c2;
  }
  return 0;
}

#endif
