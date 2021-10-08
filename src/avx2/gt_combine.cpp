#include <fiuncho/GenotypeTable.h>
#include <immintrin.h>

template <>
void GenotypeTable<uint64_t>::combine(const GenotypeTable<uint64_t> &t1,
                                      const GenotypeTable<uint64_t> &t2,
                                      GenotypeTable<uint64_t> &out) noexcept
{
    size_t i, j, k;
    // Compute bit tables for cases
    for (i = 0; i < t1.size; i += 3) {
        for (j = 0; j < 3; j++) {
            for (k = 0; k < t1.cases_words; k += 4) {
                __m256i y0 = _mm256_load_si256(
                    (__m256i *)(t2.cases + j * t1.cases_words + k));
                __m256i y1 = _mm256_load_si256(
                    (__m256i *)(t1.cases + (i + 0) * t1.cases_words + k));
                __m256i y2 = _mm256_load_si256(
                    (__m256i *)(t1.cases + (i + 1) * t1.cases_words + k));
                __m256i y3 = _mm256_load_si256(
                    (__m256i *)(t1.cases + (i + 2) * t1.cases_words + k));
                __m256i y4 = _mm256_and_si256(y0, y1);
                __m256i y5 = _mm256_and_si256(y0, y2);
                __m256i y6 = _mm256_and_si256(y0, y3);
                _mm256_store_si256(
                    (__m256i *)(out.cases + ((i + 0) * 3 + j) * t1.cases_words +
                                k),
                    y4);
                _mm256_store_si256(
                    (__m256i *)(out.cases + ((i + 1) * 3 + j) * t1.cases_words +
                                k),
                    y5);
                _mm256_store_si256(
                    (__m256i *)(out.cases + ((i + 2) * 3 + j) * t1.cases_words +
                                k),
                    y6);
            }
        }
    }
    // Compute bit tables for ctrls
    for (i = 0; i < t1.size; i += 3) {
        for (j = 0; j < 3; j++) {
            for (k = 0; k < t1.ctrls_words; k += 4) {
                __m256i y0 = _mm256_load_si256(
                    (__m256i *)(t2.ctrls + j * t1.ctrls_words + k));
                __m256i y1 = _mm256_load_si256(
                    (__m256i *)(t1.ctrls + (i + 0) * t1.ctrls_words + k));
                __m256i y2 = _mm256_load_si256(
                    (__m256i *)(t1.ctrls + (i + 1) * t1.ctrls_words + k));
                __m256i y3 = _mm256_load_si256(
                    (__m256i *)(t1.ctrls + (i + 2) * t1.ctrls_words + k));
                __m256i y4 = _mm256_and_si256(y0, y1);
                __m256i y5 = _mm256_and_si256(y0, y2);
                __m256i y6 = _mm256_and_si256(y0, y3);
                _mm256_store_si256(
                    (__m256i *)(out.ctrls + ((i + 0) * 3 + j) * t1.ctrls_words +
                                k),
                    y4);
                _mm256_store_si256(
                    (__m256i *)(out.ctrls + ((i + 1) * 3 + j) * t1.ctrls_words +
                                k),
                    y5);
                _mm256_store_si256(
                    (__m256i *)(out.ctrls + ((i + 2) * 3 + j) * t1.ctrls_words +
                                k),
                    y6);
            }
        }
    }
}
