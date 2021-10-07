#include <fiuncho/GenotypeTable.h>
#include <immintrin.h>

inline void popcnt_aux(const uint64_t *gt_tbl1, const size_t gt_size,
                       const uint64_t *gt_tbl2, const size_t words,
                       uint32_t *ct_tbl, const size_t ct_size)
{
    size_t i, j, k;
    for (i = 0; i < gt_size; ++i) {
        for (j = 0; j < 3; ++j) {
            uint32_t cnt[4] = {0, 0, 0, 0};
            for (k = 0; k < words; k += 4) {
                __m256i w1 =
                    _mm256_load_si256((__m256i *)(gt_tbl1 + i * words + k));
                __m256i w2 =
                    _mm256_load_si256((__m256i *)(gt_tbl2 + j * words + k));
                __m256i inters = _mm256_and_si256(w1, w2);
                cnt[0] += _mm_popcnt_u64(inters[0]);
                cnt[1] += _mm_popcnt_u64(inters[1]);
                cnt[2] += _mm_popcnt_u64(inters[2]);
                cnt[3] += _mm_popcnt_u64(inters[3]);
            }
            ct_tbl[i * 3 + j] = cnt[0] + cnt[1] + cnt[2] + cnt[3];
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
    popcnt_aux(t1.ctrls, t1.size, t2.ctrls, t2.ctrls_words, out.ctrls,
               out.size);
}
