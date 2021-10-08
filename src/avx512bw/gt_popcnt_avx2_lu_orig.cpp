#include <fiuncho/GenotypeTable.h>
#include <immintrin.h>

inline void popcnt_aux(const uint64_t *gt_tbl1, const size_t gt_size,
                       const uint64_t *gt_tbl2, const size_t words,
                       uint32_t *ct_tbl, const size_t ct_size)
{
    const __m256i lookup = _mm256_setr_epi8(
        /* 0 */ 0, /* 1 */ 1, /* 2 */ 1, /* 3 */ 2,
        /* 4 */ 1, /* 5 */ 2, /* 6 */ 2, /* 7 */ 3,
        /* 8 */ 1, /* 9 */ 2, /* a */ 2, /* b */ 3,
        /* c */ 2, /* d */ 3, /* e */ 3, /* f */ 4,

        /* 0 */ 0, /* 1 */ 1, /* 2 */ 1, /* 3 */ 2,
        /* 4 */ 1, /* 5 */ 2, /* 6 */ 2, /* 7 */ 3,
        /* 8 */ 1, /* 9 */ 2, /* a */ 2, /* b */ 3,
        /* c */ 2, /* d */ 3, /* e */ 3, /* f */ 4);

    const __m256i low_mask = _mm256_set1_epi8(0x0f);

    size_t i, j, k, l;
    for (i = 0; i < gt_size; ++i) {
        for (j = 0; j < 3; ++j) {
            __m256i acc = _mm256_setzero_si256();
            k = 0;
            while (k < words) {
                __m256i local = _mm256_setzero_si256();
                for (l = 0; l < 255 / 8 && k < words; l += 2, k += 8) {
                    const __m512i w1 =
                        _mm512_load_si512((__m512i *)(gt_tbl1 + i * words + k));
                    const __m512i w2 =
                        _mm512_load_si512((__m512i *)(gt_tbl2 + j * words + k));
                    const __m512i inter = _mm512_and_si512(w1, w2);
                    const __m256i lower = _mm512_extracti64x4_epi64(inter, 0);
                    __m256i lo = _mm256_and_si256(lower, low_mask);
                    __m256i hi =
                        _mm256_and_si256(_mm256_srli_epi16(lower, 4), low_mask);
                    __m256i popcnt1 = _mm256_shuffle_epi8(lookup, lo);
                    __m256i popcnt2 = _mm256_shuffle_epi8(lookup, hi);
                    const __m256i sum1 = _mm256_add_epi8(popcnt1, popcnt2);
                    const __m256i upper = _mm512_extracti64x4_epi64(inter, 1);
                    lo = _mm256_and_si256(upper, low_mask);
                    hi =
                        _mm256_and_si256(_mm256_srli_epi16(upper, 4), low_mask);
                    popcnt1 = _mm256_shuffle_epi8(lookup, lo);
                    popcnt2 = _mm256_shuffle_epi8(lookup, hi);
                    const __m256i sum2 = _mm256_add_epi8(popcnt1, popcnt2);
                    local = _mm256_add_epi8(local, sum1);
                    local = _mm256_add_epi8(local, sum2);
                }
                acc = _mm256_add_epi64(
                    acc, _mm256_sad_epu8(local, _mm256_setzero_si256()));
            }
            ct_tbl[i * 3 + j] =
                static_cast<uint64_t>(_mm256_extract_epi64(acc, 0)) +
                static_cast<uint64_t>(_mm256_extract_epi64(acc, 1)) +
                static_cast<uint64_t>(_mm256_extract_epi64(acc, 2)) +
                static_cast<uint64_t>(_mm256_extract_epi64(acc, 3));
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
