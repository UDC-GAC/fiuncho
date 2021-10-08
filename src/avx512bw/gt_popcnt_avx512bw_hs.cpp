#include <fiuncho/GenotypeTable.h>
#include <immintrin.h>

__m512i popcount(const __m512i v)
{
    const __m512i m1 = _mm512_set1_epi8(0x55);
    const __m512i m2 = _mm512_set1_epi8(0x33);
    const __m512i m4 = _mm512_set1_epi8(0x0F);

    const __m512i t1 = _mm512_sub_epi8(v, (_mm512_srli_epi16(v, 1) & m1));
    const __m512i t2 =
        _mm512_add_epi8(t1 & m2, (_mm512_srli_epi16(t1, 2) & m2));
    const __m512i t3 = _mm512_add_epi8(t2, _mm512_srli_epi16(t2, 4)) & m4;
    return _mm512_sad_epu8(t3, _mm512_setzero_si512());
}

void CSA(__m512i &h, __m512i &l, __m512i a, __m512i b, __m512i c)
{
    /*
        c b a | l h
        ------+----
        0 0 0 | 0 0
        0 0 1 | 1 0
        0 1 0 | 1 0
        0 1 1 | 0 1
        1 0 0 | 1 0
        1 0 1 | 0 1
        1 1 0 | 0 1
        1 1 1 | 1 1
        l - digit
        h - carry
    */

    l = _mm512_ternarylogic_epi32(c, b, a, 0x96);
    h = _mm512_ternarylogic_epi32(c, b, a, 0xe8);
}

inline void popcnt_harley_seal(const uint64_t *gt_tbl1, const size_t gt_size,
                               const uint64_t *gt_tbl2, const size_t words,
                               uint32_t *ct_tbl, const size_t ct_size)
{
    size_t i, j, k, l;
    // Combine two genotype tables and save its contingency table
    for (i = 0; i < gt_size; ++i) {
        for (j = 0; j < 3; ++j) {
            __m512i total = _mm512_setzero_si512();
            __m512i ones = _mm512_setzero_si512();
            __m512i twos = _mm512_setzero_si512();
            __m512i fours = _mm512_setzero_si512();
            __m512i eights = _mm512_setzero_si512();
            __m512i sixteens = _mm512_setzero_si512();
            __m512i twosA, twosB, foursA, foursB, eightsA, eightsB;

            const size_t limit = words - words % 128;
            for (k = 0; k < limit; k += 128) {
                __m512i w11 =
                    _mm512_load_si512((__m512i *)(gt_tbl1 + i * words + k + 0));
                __m512i w21 =
                    _mm512_load_si512((__m512i *)(gt_tbl2 + j * words + k + 0));
                __m512i i1 = _mm512_and_si512(w11, w21);
                __m512i w12 =
                    _mm512_load_si512((__m512i *)(gt_tbl1 + i * words + k + 8));
                __m512i w22 =
                    _mm512_load_si512((__m512i *)(gt_tbl2 + j * words + k + 8));
                __m512i i2 = _mm512_and_si512(w12, w22);
                CSA(twosA, ones, ones, i1, i2);
                w11 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + k + 16));
                w21 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + k + 16));
                i1 = _mm512_and_si512(w11, w21);
                w12 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + k + 24));
                w22 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + k + 24));
                i2 = _mm512_and_si512(w12, w22);
                CSA(twosB, ones, ones, i1, i2);
                CSA(foursA, twos, twos, twosA, twosB);
                w11 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + k + 32));
                w21 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + k + 32));
                i1 = _mm512_and_si512(w11, w21);
                w12 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + k + 40));
                w22 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + k + 40));
                i2 = _mm512_and_si512(w12, w22);
                CSA(twosA, ones, ones, i1, i2);
                w11 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + k + 48));
                w21 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + k + 48));
                i1 = _mm512_and_si512(w11, w21);
                w12 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + k + 56));
                w22 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + k + 56));
                i2 = _mm512_and_si512(w12, w22);
                CSA(twosB, ones, ones, i1, i2);
                CSA(foursB, twos, twos, twosA, twosB);
                CSA(eightsA, fours, fours, foursA, foursB);
                w11 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + k + 64));
                w21 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + k + 64));
                i1 = _mm512_and_si512(w11, w21);
                w12 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + k + 72));
                w22 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + k + 72));
                i2 = _mm512_and_si512(w12, w22);
                CSA(twosA, ones, ones, i1, i2);
                w11 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + k + 80));
                w21 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + k + 80));
                i1 = _mm512_and_si512(w11, w21);
                w12 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + k + 88));
                w22 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + k + 88));
                i2 = _mm512_and_si512(w12, w22);
                CSA(twosB, ones, ones, i1, i2);
                CSA(foursA, twos, twos, twosA, twosB);
                w11 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + k + 96));
                w21 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + k + 96));
                i1 = _mm512_and_si512(w11, w21);
                w12 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + k + 104));
                w22 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + k + 104));
                i2 = _mm512_and_si512(w12, w22);
                CSA(twosA, ones, ones, i1, i2);
                w11 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + k + 112));
                w21 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + k + 112));
                i1 = _mm512_and_si512(w11, w21);
                w12 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + k + 120));
                w22 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + k + 120));
                i2 = _mm512_and_si512(w12, w22);
                CSA(twosB, ones, ones, i1, i2);
                CSA(foursB, twos, twos, twosA, twosB);
                CSA(eightsB, fours, fours, foursA, foursB);
                CSA(sixteens, eights, eights, eightsA, eightsB);

                total = _mm512_add_epi64(total, popcount(sixteens));
            }

            total = _mm512_slli_epi64(total, 4); // * 16
            total = _mm512_add_epi64(
                total, _mm512_slli_epi64(popcount(eights), 3)); // += 8 * ...
            total = _mm512_add_epi64(
                total, _mm512_slli_epi64(popcount(fours), 2)); // += 4 * ...
            total = _mm512_add_epi64(
                total, _mm512_slli_epi64(popcount(twos), 1)); // += 2 * ...
            total = _mm512_add_epi64(total, popcount(ones));

            for (; k < words; k += 8) {
                __m512i w1 =
                    _mm512_load_si512((__m512i *)(gt_tbl1 + i * words + k));
                __m512i w2 =
                    _mm512_load_si512((__m512i *)(gt_tbl2 + j * words + k));
                __m512i inters = _mm512_and_si512(w1, w2);
                total = _mm512_add_epi64(total, popcount(inters));
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
    popcnt_harley_seal(t1.cases, t1.size, t2.cases, t1.cases_words, out.cases,
                       out.size);
    popcnt_harley_seal(t1.ctrls, t1.size, t2.ctrls, t1.ctrls_words, out.ctrls,
                       out.size);
}
