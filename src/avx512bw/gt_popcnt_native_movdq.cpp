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
            for (k = 0; k < words; k += 8) {
                __m512i w1 =
                    _mm512_load_si512((__m512i *)(gt_tbl1 + i * words + k));
                __m512i w2 =
                    _mm512_load_si512((__m512i *)(gt_tbl2 + j * words + k));
                __m512i inters = _mm512_and_si512(w1, w2);
                __m256i lower256 = _mm512_extracti64x4_epi64(inters, 0);
                __m128i lower128 = _mm256_extracti64x2_epi64(lower256, 0);
                uint64_t lower64 = _mm_cvtsi128_si64(lower128);
                cnt += _mm_popcnt_u64(lower64);
                __m128i temp2 =
                    (__m128i)_mm_movehl_ps((__m128)lower128, (__m128)lower128);
                uint64_t upper64 = _mm_cvtsi128_si64(temp2);
                cnt += _mm_popcnt_u64(upper64);
                __m128i upper128 = _mm256_extracti64x2_epi64(lower256, 1);
                lower64 = _mm_cvtsi128_si64(upper128);
                cnt += _mm_popcnt_u64(lower64);
                temp2 =
                    (__m128i)_mm_movehl_ps((__m128)upper128, (__m128)upper128);
                upper64 = _mm_cvtsi128_si64(temp2);
                cnt += _mm_popcnt_u64(upper64);
                __m256i upper256 = _mm512_extracti64x4_epi64(inters, 1);
                lower128 = _mm256_extracti64x2_epi64(upper256, 0);
                lower64 = _mm_cvtsi128_si64(lower128);
                cnt += _mm_popcnt_u64(lower64);
                temp2 =
                    (__m128i)_mm_movehl_ps((__m128)lower128, (__m128)lower128);
                upper64 = _mm_cvtsi128_si64(temp2);
                cnt += _mm_popcnt_u64(upper64);
                upper128 = _mm256_extracti64x2_epi64(upper256, 1);
                lower64 = _mm_cvtsi128_si64(upper128);
                cnt += _mm_popcnt_u64(lower64);
                temp2 =
                    (__m128i)_mm_movehl_ps((__m128)upper128, (__m128)upper128);
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
