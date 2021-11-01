#include <fiuncho/GenotypeTable.h>
#include <immintrin.h>

inline void combine_aux(const uint64_t *gt_tbl1, const size_t size1,
                        const uint64_t *gt_tbl2, const size_t words,
                        uint64_t *gt_tbl3)
{
    size_t i, j, k;
    // Combine two genotype tables and store it in a larger table
    // Inner gt_tbl2 loop is unrolled by its size
    const __m256i *ptr1 = (const __m256i *)gt_tbl1;
    for (i = 0; i < size1; ++i) {
        const __m256i *ptr2_1 = (const __m256i *)(gt_tbl2 + 0 * words);
        const __m256i *ptr2_2 = (const __m256i *)(gt_tbl2 + 1 * words);
        const __m256i *ptr2_3 = (const __m256i *)(gt_tbl2 + 2 * words);
        __m256i *ptr3_1 = (__m256i *)(gt_tbl3 + (i * 3 + 0) * words);
        __m256i *ptr3_2 = (__m256i *)(gt_tbl3 + (i * 3 + 1) * words);
        __m256i *ptr3_3 = (__m256i *)(gt_tbl3 + (i * 3 + 2) * words);
        for (k = 0; k < words; k += 4) {
            __m256i y0 = _mm256_load_si256(ptr1++);
            __m256i y1 = _mm256_load_si256(ptr2_1++);
            __m256i y2 = _mm256_load_si256(ptr2_2++);
            __m256i y3 = _mm256_load_si256(ptr2_3++);
            _mm256_store_si256(ptr3_1++, _mm256_and_si256(y0, y1));
            _mm256_store_si256(ptr3_2++, _mm256_and_si256(y0, y2));
            _mm256_store_si256(ptr3_3++, _mm256_and_si256(y0, y3));
        }
    }
}

template <>
void GenotypeTable<uint64_t>::combine(const GenotypeTable<uint64_t> &t1,
                                      const GenotypeTable<uint64_t> &t2,
                                      GenotypeTable<uint64_t> &out) noexcept
{
    combine_aux(t1.cases, t1.size, t2.cases, t1.cases_words, out.cases);
    combine_aux(t1.ctrls, t1.size, t2.ctrls, t1.ctrls_words, out.ctrls);
}
