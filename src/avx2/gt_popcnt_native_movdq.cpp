#include <fiuncho/GenotypeTable.h>
#include <immintrin.h>

inline void popcnt_aux(const uint64_t *gt_tbl1, const size_t gt_size,
                       const uint64_t *gt_tbl2, const size_t words,
                       uint32_t *ct_tbl, const size_t ct_size)
{
    size_t i, j, k;
    for (i = 0; i < gt_size; ++i) {
        for (j = 0; j < 3; ++j) {
            uint32_t cnt = 0;
            for (k = 0; k < words; k += 4) {
                __m256i w1 =
                    _mm256_load_si256((__m256i *)(gt_tbl1 + i * words + k));
                __m256i w2 =
                    _mm256_load_si256((__m256i *)(gt_tbl2 + j * words + k));
                __m256i inters = _mm256_and_si256(w1, w2);
                const __m128i lo = _mm256_extracti64x2_epi64(inters, 0);
                uint64_t lower64 = _mm_cvtsi128_si64(lo);
                cnt += _mm_popcnt_u64(lower64);
                __m128i temp2 = (__m128i)_mm_movehl_ps((__m128)lo, (__m128)lo);
                uint64_t upper64 = _mm_cvtsi128_si64(temp2);
                cnt += _mm_popcnt_u64(upper64);
                const __m128i hi = _mm256_extracti64x2_epi64(inters, 1);
                lower64 = _mm_cvtsi128_si64(hi);
                cnt += _mm_popcnt_u64(lower64);
                temp2 = (__m128i)_mm_movehl_ps((__m128)hi, (__m128)hi);
                upper64 = _mm_cvtsi128_si64(temp2);
                cnt += _mm_popcnt_u64(upper64);
            }
            ct_tbl[i * 3 + j] = cnt;
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
