#include <fiuncho/GenotypeTable.h>
#include <immintrin.h>

inline void combine_aux(const uint64_t *gt_tbl1, const size_t gt_size,
                        const uint64_t *gt_tbl2, const size_t words,
                        uint64_t *gt_tbl3)
{
    size_t i, j, k;
    // Combine two genotype tables and store it in a larger table
    // Outer loop is unrolled by a factor of 3
    for (i = 0; i < gt_size; i += 3) {
        for (j = 0; j < 3; j++) {
            for (k = 0; k < words; k += 8) {
                __m512i z0 =
                    _mm512_load_si512((__m512i *)(gt_tbl2 + j * words + k));
                __m512i z1 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + (i + 0) * words + k));
                __m512i z2 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + (i + 1) * words + k));
                __m512i z3 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + (i + 2) * words + k));
                __m512i z4 = _mm512_and_si512(z0, z1);
                __m512i z5 = _mm512_and_si512(z0, z2);
                __m512i z6 = _mm512_and_si512(z0, z3);
                _mm512_store_si512(
                    (__m512i *)(gt_tbl3 + ((i + 0) * 3 + j) * words + k), z4);
                _mm512_store_si512(
                    (__m512i *)(gt_tbl3 + ((i + 1) * 3 + j) * words + k), z5);
                _mm512_store_si512(
                    (__m512i *)(gt_tbl3 + ((i + 2) * 3 + j) * words + k), z6);
            }
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
