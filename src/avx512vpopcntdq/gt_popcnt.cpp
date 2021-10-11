#include <fiuncho/GenotypeTable.h>
#include <immintrin.h>

inline void popcnt(const uint64_t *gt_tbl1, const size_t gt_size,
                   const uint64_t *gt_tbl2, const size_t words,
                   uint32_t *ct_tbl, const size_t ct_size)
{
    size_t i, j, k, l;
    // Combine two genotype tables and save its contingency table
    for (i = 0; i < gt_size; ++i) {
        for (j = 0; j < 3; ++j) {
            __m512i total = _mm512_setzero_si512();
            for (k = 0; k < words; k += 8) {
                __m512i w1 =
                    _mm512_load_si512((__m512i *)(gt_tbl1 + i * words + k));
                __m512i w2 =
                    _mm512_load_si512((__m512i *)(gt_tbl2 + j * words + k));
                __m512i inters = _mm512_and_si512(w1, w2);
                total = _mm512_add_epi64(total, _mm512_popcnt_epi64(inters));
            }
            const __m256i lo = _mm512_extracti64x4_epi64(total, 0);
            const __m256i hi = _mm512_extracti64x4_epi64(total, 1);
            ct_tbl[i * 3 + j] =
                static_cast<uint64_t>(_mm256_extract_epi64(lo, 0)) +
                static_cast<uint64_t>(_mm256_extract_epi64(lo, 1)) +
                static_cast<uint64_t>(_mm256_extract_epi64(lo, 2)) +
                static_cast<uint64_t>(_mm256_extract_epi64(lo, 3)) +
                static_cast<uint64_t>(_mm256_extract_epi64(hi, 0)) +
                static_cast<uint64_t>(_mm256_extract_epi64(hi, 1)) +
                static_cast<uint64_t>(_mm256_extract_epi64(hi, 2)) +
                static_cast<uint64_t>(_mm256_extract_epi64(hi, 3));
        }
    }
    for (i = gt_size * 3; i < ct_size; ++i) {
        ct_tbl[i] = 0;
    }
}

template <>
template <>
void GenotypeTable<uint64_t>::combine_and_popcnt(
    const GenotypeTable<uint64_t> &t1, const GenotypeTable<uint64_t> &t2,
    ContingencyTable<uint32_t> &out) noexcept
{
    popcnt(t1.cases, t1.size, t2.cases, t1.cases_words, out.cases, out.size);
    popcnt(t1.ctrls, t1.size, t2.ctrls, t1.ctrls_words, out.ctrls, out.size);
}
