#include <fiuncho/GenotypeTable.h>
#include <immintrin.h>

inline void iter(const uint64_t *ptr1, const uint64_t *ptr2,
                 const __m256i &lookup, const __m256i &low_mask, __m256i &local)
{
    __m256i o1 = _mm256_load_si256((__m256i *)ptr1);
    __m256i o2 = _mm256_load_si256((__m256i *)ptr2);
    const __m256i vec = _mm256_and_si256(o1, o2);
    const __m256i lo = _mm256_and_si256(vec, low_mask);
    const __m256i hi = _mm256_and_si256(_mm256_srli_epi16(vec, 4), low_mask);
    const __m256i popcnt1 = _mm256_shuffle_epi8(lookup, lo);
    const __m256i popcnt2 = _mm256_shuffle_epi8(lookup, hi);
    local = _mm256_add_epi8(local, popcnt1);
    local = _mm256_add_epi8(local, popcnt2);
}

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

    size_t i, j, k;
    for (i = 0; i < gt_size; ++i) {
        for (j = 0; j < 3; ++j) {
            __m256i acc = _mm256_setzero_si256();
            for (k = 0; k + 32 <= words; k += 32) {
                __m256i local = _mm256_setzero_si256();
                iter(gt_tbl1 + i * words + k + 0, gt_tbl2 + j * words + k + 0,
                     lookup, low_mask, local);
                iter(gt_tbl1 + i * words + k + 4, gt_tbl2 + j * words + k + 4,
                     lookup, low_mask, local);
                iter(gt_tbl1 + i * words + k + 8, gt_tbl2 + j * words + k + 8,
                     lookup, low_mask, local);
                iter(gt_tbl1 + i * words + k + 12, gt_tbl2 + j * words + k + 12,
                     lookup, low_mask, local);
                iter(gt_tbl1 + i * words + k + 16, gt_tbl2 + j * words + k + 16,
                     lookup, low_mask, local);
                iter(gt_tbl1 + i * words + k + 20, gt_tbl2 + j * words + k + 20,
                     lookup, low_mask, local);
                iter(gt_tbl1 + i * words + k + 24, gt_tbl2 + j * words + k + 24,
                     lookup, low_mask, local);
                iter(gt_tbl1 + i * words + k + 28, gt_tbl2 + j * words + k + 28,
                     lookup, low_mask, local);
                acc = _mm256_add_epi64(
                    acc, _mm256_sad_epu8(local, _mm256_setzero_si256()));
            }

            __m256i local = _mm256_setzero_si256();
            for (; k < words; k += 4) {
                iter(gt_tbl1 + i * words + k, gt_tbl2 + j * words + k, lookup,
                     low_mask, local);
            }
            acc = _mm256_add_epi64(
                acc, _mm256_sad_epu8(local, _mm256_setzero_si256()));

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
