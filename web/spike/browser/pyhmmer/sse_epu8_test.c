#include <stdio.h>
#include <stdint.h>
#include <smmintrin.h>
static void dump8u(const char *name, __m128i v) {
  uint8_t a[16]; _mm_storeu_si128((__m128i *)a, v);
  printf("%s:", name); for (int i = 0; i < 16; ++i) printf(" %u", a[i]); printf("\n");
}
int main(void) {
  uint8_t ba[16], bb[16];
  for (int i = 0; i < 16; ++i) { ba[i] = (uint8_t)(i * 37 + 11); bb[i] = (uint8_t)(255 - i * 13); }
  __m128i p = _mm_loadu_si128((const __m128i *)ba), q = _mm_loadu_si128((const __m128i *)bb);
  dump8u("adds_epu8", _mm_adds_epu8(p, q));
  dump8u("subs_epu8", _mm_subs_epu8(p, q));
  dump8u("max_epu8", _mm_max_epu8(p, q));
  dump8u("min_epu8", _mm_min_epu8(p, q));
  dump8u("cmpeq_epi8", _mm_cmpeq_epi8(p, q));
  dump8u("set1_epi8", _mm_set1_epi8((char)-57));
  dump8u("slli_si128_1", _mm_slli_si128(p, 1));
  dump8u("srli_si128_1", _mm_srli_si128(p, 1));
  dump8u("shuffle_epi8", _mm_shuffle_epi8(p, q));
  dump8u("alignr_epi8", _mm_alignr_epi8(p, q, 3));
  dump8u("blendv_epi8", _mm_blendv_epi8(p, q, _mm_cmpgt_epi8(p, q)));
  printf("movemask: %d\n", _mm_movemask_epi8(p));
  printf("extract_epi8_5: %d\n", _mm_extract_epi8(p, 5));
  printf("insert_epi8: "); { __m128i r = _mm_insert_epi8(p, 77, 3); dump8u("", r); }
  return 0;
}
