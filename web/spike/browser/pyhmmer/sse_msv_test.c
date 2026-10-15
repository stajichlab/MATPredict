#include <stdio.h>
#include <stdint.h>
#include <smmintrin.h>
static void dump16(const char *n, __m128i v) { int16_t a[8]; _mm_storeu_si128((__m128i *)a, v); printf("%s:", n); for (int i = 0; i < 8; ++i) printf(" %d", a[i]); printf("\n"); }
static void dump8(const char *n, __m128i v) { int8_t a[16]; _mm_storeu_si128((__m128i *)a, v); printf("%s:", n); for (int i = 0; i < 16; ++i) printf(" %d", a[i]); printf("\n"); }
static void dump32(const char *n, __m128i v) { int32_t a[4]; _mm_storeu_si128((__m128i *)a, v); printf("%s:", n); for (int i = 0; i < 4; ++i) printf(" %d", a[i]); printf("\n"); }
int main(void) {
  int8_t ba[16], bb[16];
  for (int i = 0; i < 16; ++i) { ba[i] = (int8_t)(i * 37 - 90); bb[i] = (int8_t)(120 - i * 17); }
  __m128i p = _mm_loadu_si128((const __m128i *)ba), q = _mm_loadu_si128((const __m128i *)bb);
  printf("extract_epi16_0: %d\n", _mm_extract_epi16(p, 0));
  printf("extract_epi16_3: %d\n", _mm_extract_epi16(p, 3));
  printf("extract_epi16_7: %d\n", _mm_extract_epi16(p, 7));
  dump32("shuffle_epi32_1b", _mm_shuffle_epi32(p, 0x1B));
  dump32("shuffle_epi32_e4", _mm_shuffle_epi32(q, 0xE4));
  dump16("shufflelo_epi16_1b", _mm_shufflelo_epi16(p, 0x1B));
  dump16("shufflelo_epi16_39", _mm_shufflelo_epi16(q, 0x39));
  dump8("cmpgt_epi8", _mm_cmpgt_epi8(p, q));
  return 0;
}
