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

            const uint64_t limit = words - words % 128;
            k = 0;
            for (; k < limit; k += 128) {
                __m512i tbl1_word1 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + (k + 0)));
                __m512i tbl2_word1 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + (k + 0)));
                __m512i int1 = _mm512_and_si512(tbl1_word1, tbl2_word1);
                __m512i tbl1_word2 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + (k + 8)));
                __m512i tbl2_word2 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + (k + 8)));
                __m512i int2 = _mm512_and_si512(tbl1_word2, tbl2_word2);
                CSA(twosA, ones, ones, int1, int2);
                tbl1_word1 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + (k + 16)));
                tbl2_word1 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + (k + 16)));
                int1 = _mm512_and_si512(tbl1_word1, tbl2_word1);
                tbl1_word2 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + (k + 24)));
                tbl2_word2 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + (k + 24)));
                int2 = _mm512_and_si512(tbl1_word2, tbl2_word2);
                CSA(twosB, ones, ones, int1, int2);
                CSA(foursA, twos, twos, twosA, twosB);
                tbl1_word1 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + (k + 32)));
                tbl2_word1 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + (k + 32)));
                int1 = _mm512_and_si512(tbl1_word1, tbl2_word1);
                tbl1_word2 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + (k + 40)));
                tbl2_word2 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + (k + 40)));
                int2 = _mm512_and_si512(tbl1_word2, tbl2_word2);
                CSA(twosA, ones, ones, int1, int2);
                tbl1_word1 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + (k + 48)));
                tbl2_word1 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + (k + 48)));
                int1 = _mm512_and_si512(tbl1_word1, tbl2_word1);
                tbl1_word2 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + (k + 56)));
                tbl2_word2 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + (k + 56)));
                int2 = _mm512_and_si512(tbl1_word2, tbl2_word2);
                CSA(twosB, ones, ones, int1, int2);
                CSA(foursB, twos, twos, twosA, twosB);
                CSA(eightsA, fours, fours, foursA, foursB);
                tbl1_word1 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + (k + 64)));
                tbl2_word1 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + (k + 64)));
                int1 = _mm512_and_si512(tbl1_word1, tbl2_word1);
                tbl1_word2 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + (k + 72)));
                tbl2_word2 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + (k + 72)));
                int2 = _mm512_and_si512(tbl1_word2, tbl2_word2);
                CSA(twosA, ones, ones, int1, int2);
                tbl1_word1 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + (k + 80)));
                tbl2_word1 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + (k + 80)));
                int1 = _mm512_and_si512(tbl1_word1, tbl2_word1);
                tbl1_word2 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + (k + 88)));
                tbl2_word2 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + (k + 88)));
                int2 = _mm512_and_si512(tbl1_word2, tbl2_word2);
                CSA(twosB, ones, ones, int1, int2);
                CSA(foursA, twos, twos, twosA, twosB);
                tbl1_word1 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + (k + 96)));
                tbl2_word1 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + (k + 96)));
                int1 = _mm512_and_si512(tbl1_word1, tbl2_word1);
                tbl1_word2 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + (k + 104)));
                tbl2_word2 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + (k + 104)));
                int2 = _mm512_and_si512(tbl1_word2, tbl2_word2);
                CSA(twosA, ones, ones, int1, int2);
                tbl1_word1 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + (k + 112)));
                tbl2_word1 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + (k + 112)));
                int1 = _mm512_and_si512(tbl1_word1, tbl2_word1);
                tbl1_word2 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + (k + 120)));
                tbl2_word2 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + (k + 120)));
                int2 = _mm512_and_si512(tbl1_word2, tbl2_word2);
                CSA(twosB, ones, ones, int1, int2);
                CSA(foursB, twos, twos, twosA, twosB);
                CSA(eightsB, fours, fours, foursA, foursB);
                CSA(sixteens, eights, eights, eightsA, eightsB);

                total = _mm512_add_epi64(total, popcount(sixteens));
            }

            for (; k < words; k++) {
                __m512i tbl1_word =
                    _mm512_load_si512((__m512i *)(gt_tbl1 + i * words + k));
                __m512i tbl2_word =
                    _mm512_load_si512((__m512i *)(gt_tbl2 + j * words + k));
                __m512i inters = _mm512_and_si512(tbl1_word, tbl2_word);
                total = _mm512_add_epi64(total, popcount(inters));
            }

            ct_tbl[i * 3 + j] = _mm512_reduce_add_epi64(total);
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
    popcnt_harley_seal(t1.ctrls, t1.size, t2.ctrls, t2.ctrls_words, out.ctrls,
                       out.size);
}
