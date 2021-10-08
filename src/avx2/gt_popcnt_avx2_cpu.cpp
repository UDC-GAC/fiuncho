#include <fiuncho/GenotypeTable.h>
#include <immintrin.h>
#include <x86intrin.h>

inline void iter(const uint64_t *ptr1, const uint64_t *ptr2, uint32_t &result)
{
    __m256i w1 = _mm256_load_si256((__m256i *)ptr1);
    __m256i w2 = _mm256_load_si256((__m256i *)ptr2);
    __m256i inter = _mm256_and_si256(w1, w2);
    result += _popcnt64(_mm256_extract_epi64(inter, 0));
    result += _popcnt64(_mm256_extract_epi64(inter, 1));
    result += _popcnt64(_mm256_extract_epi64(inter, 2));
    result += _popcnt64(_mm256_extract_epi64(inter, 3));
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
                iter(gt_tbl1 + i * words + k + 0, gt_tbl2 + j * words + k + 0,
                     result);
                iter(gt_tbl1 + i * words + k + 4, gt_tbl2 + j * words + k + 4,
                     result);
                iter(gt_tbl1 + i * words + k + 8, gt_tbl2 + j * words + k + 8,
                     result);
                iter(gt_tbl1 + i * words + k + 12, gt_tbl2 + j * words + k + 12,
                     result);
            }
            for (; k < words; k += 4) {
                iter(gt_tbl1 + i * words + k, gt_tbl2 + j * words + k, result);
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
