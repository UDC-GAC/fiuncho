#include <cmath>
#include <fiuncho/GenotypeTable.h>
#include <immintrin.h>

template <>
GenotypeTable<uint64_t>::GenotypeTable(const short order,
                                       const size_t cases_words,
                                       const size_t ctrls_words)
    : order(order), size(std::pow(3, order)), cases_words(cases_words),
      ctrls_words(ctrls_words), alloc(std::make_unique<uint64_t[]>(
                                    size * (cases_words + ctrls_words) + 8)),
      cases((uint64_t *)((((uintptr_t)alloc.get()) + 63) / 64 * 64)),
      ctrls(cases + size * cases_words)
{
}

inline void combine_aux(const uint64_t *gt_tbl1, const size_t gt_size,
                        const uint64_t *gt_tbl2, const size_t words,
                        uint64_t *gt_tbl3)
{
    size_t i, j, k;
    // Combine two genotype tables and store it in a larger table
    // Outer loop is unrolled by a factor of 3
    for (i = 0; i < gt_size; i += 3) {
        for (j = 0; j < 3; j++) {
            for (k = 0; k < words; k += 8) {
                __m512i z0 =
                    _mm512_load_si512((__m512i *)(gt_tbl2 + j * words + k));
                __m512i z1 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + (i + 0) * words + k));
                __m512i z2 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + (i + 1) * words + k));
                __m512i z3 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + (i + 2) * words + k));
                __m512i z4 = _mm512_and_si512(z0, z1);
                __m512i z5 = _mm512_and_si512(z0, z2);
                __m512i z6 = _mm512_and_si512(z0, z3);
                _mm512_store_si512(
                    (__m512i *)(gt_tbl3 + ((i + j) * 3 + 0) * words + k), z4);
                _mm512_store_si512(
                    (__m512i *)(gt_tbl3 + ((i + j) * 3 + 1) * words + k), z5);
                _mm512_store_si512(
                    (__m512i *)(gt_tbl3 + ((i + j) * 3 + 2) * words + k), z6);
            }
        }
    }
}

template <>
void GenotypeTable<uint64_t>::combine(const GenotypeTable<uint64_t> &t1,
                                      const GenotypeTable<uint64_t> &t2,
                                      GenotypeTable<uint64_t> &out) noexcept
{
    combine_aux(t1.cases, t1.size, t2.cases, t1.cases_words, out.cases);
    combine_aux(t1.ctrls, t1.size, t2.ctrls, t1.ctrls_words, out.ctrls);
}

const __m512i lookup = _mm512_setr_epi64(
    0x0302020102010100llu, 0x0403030203020201llu, 0x0302020102010100llu,
    0x0403030203020201llu, 0x0302020102010100llu, 0x0403030203020201llu,
    0x0302020102010100llu, 0x0403030203020201llu);
const __m512i low_mask = _mm512_set1_epi8(0x0f);

inline void popcnt_aux(const uint64_t *gt_tbl1, const size_t gt_size,
                       const uint64_t *gt_tbl2, const size_t words,
                       uint32_t *ct_tbl, const size_t ct_size)
{
    size_t i, j, k, l;
    // Combine two genotype tables and save its contingency table
    for (i = 0; i < gt_size; ++i) {
        for (j = 0; j < 3; ++j) {
            k = 0;
            __m512i acc = _mm512_setzero_si512();
            while (k < words) {
                // 512-bit value holding the number of individuals counted in
                // each 8-bit value inside the 512-bit word
                __m512i local = _mm512_setzero_si512();
                // Each 8-bit integer can represent a value in [0, 255], and
                // each iteration can increase the count by 8, so the maximum
                // number of iterations before clearing the buffer is 255/8
                for (l = 0; l < 255 / 8 && k < words; ++l, k += 8) {
                    // Read values to combine from input genotype tables
                    __m512i z0 =
                        _mm512_load_si512((__m512i *)(gt_tbl2 + j * words + k));
                    __m512i z1 =
                        _mm512_load_si512((__m512i *)(gt_tbl1 + i * words + k));
                    // Find the intersection of individuals
                    __m512i z2 = _mm512_and_si512(z0, z1);
                    // Count the number of individuals on the intersection
                    const __m512i lo = _mm512_and_si512(z2, low_mask);
                    const __m512i hi =
                        _mm512_and_si512(_mm512_srli_epi32(z2, 4), low_mask);
                    const __m512i popcnt1 = _mm512_shuffle_epi8(lookup, lo);
                    const __m512i popcnt2 = _mm512_shuffle_epi8(lookup, hi);
                    local = _mm512_add_epi8(local, popcnt1);
                    local = _mm512_add_epi8(local, popcnt2);
                }
                // Once the 512-bit buffer holding the 8-bit integers cannot
                // safely count more individuals, or there are no more words to
                // process, accumulate the results in a 512-bit buffer
                // containing 64-bit integers
                acc = _mm512_add_epi64(
                    acc, _mm512_sad_epu8(local, _mm512_setzero_si512()));
            }
            // Once all the words of a row have been processed, reduce the sum
            // into a single 64-bit value and store it in the contingency table
            ct_tbl[i * 3 + j] = _mm512_reduce_add_epi64(acc);
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
