/* Differential test: the SSE intrinsics used by miniprot's nasw-sse.c, evaluated on fixed
 * inputs. Build natively with gcc -msse4.1 (reference) and with emcc -msse4.1 -msimd128
 * (WASM). The outputs must match exactly. */
#include <stdio.h>
#include <stdint.h>
#include <smmintrin.h>

static void dump16(const char *name, __m128i v) {
	int16_t a[8];
	_mm_storeu_si128((__m128i *)a, v);
	printf("%s:", name);
	for (int i = 0; i < 8; ++i) printf(" %d", a[i]);
	printf("\n");
}

static void dump32(const char *name, __m128i v) {
	int32_t a[4];
	_mm_storeu_si128((__m128i *)a, v);
	printf("%s:", name);
	for (int i = 0; i < 4; ++i) printf(" %d", a[i]);
	printf("\n");
}

static void dump8(const char *name, __m128i v) {
	uint8_t a[16];
	_mm_storeu_si128((__m128i *)a, v);
	printf("%s:", name);
	for (int i = 0; i < 16; ++i) printf(" %u", a[i]);
	printf("\n");
}

int main(void) {
	/* Fixed, mixed-sign inputs, including saturation edges. */
	int16_t ia[8] = {  -32768, 32767, -1, 0, 1234, -4321, 7, -7 };
	int16_t ib[8] = {   5, -5, 32767, -32768, 1, 1, -100, 100 };
	int32_t ja[4] = { -2147483647 - 1, 2147483647, -9, 42 };
	int32_t jb[4] = { 3, -3, 2147483647, -2147483647 - 1 };
	uint8_t ba[16], bb[16];
	for (int i = 0; i < 16; ++i) { ba[i] = (uint8_t)(i * 37 + 11); bb[i] = (uint8_t)(255 - i * 13); }

	__m128i a = _mm_loadu_si128((const __m128i *)ia);
	__m128i b = _mm_loadu_si128((const __m128i *)ib);
	__m128i x = _mm_loadu_si128((const __m128i *)ja);
	__m128i y = _mm_loadu_si128((const __m128i *)jb);
	__m128i p = _mm_loadu_si128((const __m128i *)ba);
	__m128i q = _mm_loadu_si128((const __m128i *)bb);

	dump16("subs_epi16", _mm_subs_epi16(a, b));
	dump16("adds_epi16", _mm_adds_epi16(a, b));
	dump16("max_epi16", _mm_max_epi16(a, b));
	dump16("cmpgt_epi16", _mm_cmpgt_epi16(a, b));
	dump16("set1_epi16", _mm_set1_epi16(-321));
	dump16("insert_epi16_0", _mm_insert_epi16(a, 1234, 0));
	dump16("insert_epi16_5", _mm_insert_epi16(a, -2, 5));
	dump16("slli_si128_3", _mm_slli_si128(a, 3));
	dump16("srli_si128_5", _mm_srli_si128(a, 5));
	printf("extract_epi16_2: %d\n", _mm_extract_epi16(a, 2));
	dump16("cmpgt_epi16_cmp", _mm_cmpgt_epi16(b, a));

	dump32("max_epi32", _mm_max_epi32(x, y));
	dump32("sub_epi32", _mm_sub_epi32(x, y));
	dump32("add_epi32", _mm_add_epi32(x, y));
	dump32("cmpgt_epi32", _mm_cmpgt_epi32(x, y));
	dump32("set1_epi32", _mm_set1_epi32(-77));
	dump32("insert_epi32_2", _mm_insert_epi32(x, -5, 2));
	dump32("insert_epi32_3", _mm_insert_epi32(x, 99, 3));

	dump8("blendv_epi8", _mm_blendv_epi8(p, q, _mm_cmpgt_epi8(p, q)));
	dump8("and", _mm_and_si128(p, q));
	dump8("or", _mm_or_si128(p, q));
	dump8("andnot", _mm_andnot_si128(p, q));
	printf("movemask_epi8: %d\n", _mm_movemask_epi8(p));
	printf("movemask_cmp: %d\n", _mm_movemask_epi8(_mm_cmpgt_epi8(p, q)));
	return 0;
}
