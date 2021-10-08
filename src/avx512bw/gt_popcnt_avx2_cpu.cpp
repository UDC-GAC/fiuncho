#include <fiuncho/GenotypeTable.h>
#include <immintrin.h>
#include <x86intrin.h>

inline void iter(const __m256i &val, uint32_t &result)
{
    result += _popcnt64(_mm256_extract_epi64(val, 0));
    result += _popcnt64(_mm256_extract_epi64(val, 1));
    result += _popcnt64(_mm256_extract_epi64(val, 2));
    result += _popcnt64(_mm256_extract_epi64(val, 3));
}

inline void popcnt_aux(const uint64_t *gt_tbl1, const size_t gt_size,
                       const uint64_t *gt_tbl2, const size_t words,
                       uint32_t *ct_tbl, const size_t ct_size)
{
    size_t i, j, k;
    for (i = 0; i < gt_size; ++i) {
        for (j = 0; j < 3; ++j) {
            uint32_t result = 0;
            for (k = 0; k + 16 <= words; k += 16) {
                __m512i w1 =
                    _mm512_load_si512((__m512i *)(gt_tbl1 + i * words + k + 0));
                __m512i w2 =
                    _mm512_load_si512((__m512i *)(gt_tbl2 + j * words + k + 0));
                __m512i inter = _mm512_and_si512(w1, w2);
                __m256i lower = _mm512_extracti64x4_epi64(inter, 0);
                __m256i upper = _mm512_extracti64x4_epi64(inter, 1);
                iter(lower, result);
                iter(upper, result);
                w1 =
                    _mm512_load_si512((__m512i *)(gt_tbl1 + i * words + k + 8));
                w2 =
                    _mm512_load_si512((__m512i *)(gt_tbl2 + j * words + k + 8));
                inter = _mm512_and_si512(w1, w2);
                lower = _mm512_extracti64x4_epi64(inter, 0);
                upper = _mm512_extracti64x4_epi64(inter, 1);
                iter(lower, result);
                iter(upper, result);
            }
            for (; k < words; k += 8) {
                __m512i w1 =
                    _mm512_load_si512((__m512i *)(gt_tbl1 + i * words + k));
                __m512i w2 =
                    _mm512_load_si512((__m512i *)(gt_tbl2 + j * words + k));
                __m512i inter = _mm512_and_si512(w1, w2);
                __m256i lower = _mm512_extracti64x4_epi64(inter, 0);
                __m256i upper = _mm512_extracti64x4_epi64(inter, 1);
                iter(lower, result);
                iter(upper, result);
            }
            ct_tbl[i * 3 + j] = result;
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
    popcnt_aux(t1.cases, t1.size, t2.cases, t1.cases_words, out.cases,
               out.size);
    popcnt_aux(t1.ctrls, t1.size, t2.ctrls, t1.ctrls_words, out.ctrls,
               out.size);
}
