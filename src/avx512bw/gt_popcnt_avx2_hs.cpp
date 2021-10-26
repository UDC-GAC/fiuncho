/*
 * Copyright (c) 2008-2016, Wojciech Muła
 * Copyright (c) 2016, Kim Walisch
 * Copyright (c) 2016, Dan Luu
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <fiuncho/GenotypeTable.h>
#include <immintrin.h>

__m256i popcount(const __m256i v)
{
    const __m256i m1 = _mm256_set1_epi8(0x55);
    const __m256i m2 = _mm256_set1_epi8(0x33);
    const __m256i m4 = _mm256_set1_epi8(0x0F);

    const __m256i t1 = _mm256_sub_epi8(v, (_mm256_srli_epi16(v, 1) & m1));
    const __m256i t2 =
        _mm256_add_epi8(t1 & m2, (_mm256_srli_epi16(t1, 2) & m2));
    const __m256i t3 = _mm256_add_epi8(t2, _mm256_srli_epi16(t2, 4)) & m4;
    return _mm256_sad_epu8(t3, _mm256_setzero_si256());
}

void CSA(__m256i &h, __m256i &l, __m256i a, __m256i b, __m256i c)
{
    const __m256i u = a ^ b;
    h = (a & b) | (u & c);
    l = u ^ c;
}

/**
 * @brief Implements the combination of two genotype subtables and subsequent
 * contingency table computation. The function is a fork of the function
 * ``popcnt_AVX2_harley_seal`` from
 * https://github.com/WojciechMula/sse-popcount/tree/6feb3dba32c526b17de01e931c116900e3a23104,
 * modified to interleave the genotype table combination operations with the
 * `popcount` operations.
 */

inline void
avx512bw_popcnt_256_harley_seal(const uint64_t *gt_tbl1, const size_t gt_size,
                                const uint64_t *gt_tbl2, const size_t words,
                                uint32_t *ct_tbl, const size_t ct_size)
{
    size_t i, j, k, l;
    for (i = 0; i < gt_size; ++i) {
        for (j = 0; j < 3; ++j) {
            __m256i total = _mm256_setzero_si256();
            __m256i ones = _mm256_setzero_si256();
            __m256i twos = _mm256_setzero_si256();
            __m256i fours = _mm256_setzero_si256();
            __m256i eights = _mm256_setzero_si256();
            __m256i sixteens = _mm256_setzero_si256();
            __m256i twosA, twosB, foursA, foursB, eightsA, eightsB;

            const size_t limit = words - words % 64;
            for (k = 0; k < limit; k += 64) {
                __m512i w1 =
                    _mm512_load_si512((__m512i *)(gt_tbl1 + i * words + k + 0));
                __m512i w2 =
                    _mm512_load_si512((__m512i *)(gt_tbl2 + j * words + k + 0));
                __m512i inter = _mm512_and_si512(w1, w2);
                __m256i lower = _mm512_extracti64x4_epi64(inter, 0);
                __m256i upper = _mm512_extracti64x4_epi64(inter, 1);
                CSA(twosA, ones, ones, lower, upper);
                w1 =
                    _mm512_load_si512((__m512i *)(gt_tbl1 + i * words + k + 8));
                w2 =
                    _mm512_load_si512((__m512i *)(gt_tbl2 + j * words + k + 8));
                inter = _mm512_and_si512(w1, w2);
                lower = _mm512_extracti64x4_epi64(inter, 0);
                upper = _mm512_extracti64x4_epi64(inter, 1);
                CSA(twosB, ones, ones, lower, upper);
                CSA(foursA, twos, twos, twosA, twosB);
                w1 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + k + 16));
                w2 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + k + 16));
                inter = _mm512_and_si512(w1, w2);
                lower = _mm512_extracti64x4_epi64(inter, 0);
                upper = _mm512_extracti64x4_epi64(inter, 1);
                CSA(twosA, ones, ones, lower, upper);
                w1 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + k + 24));
                w2 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + k + 24));
                inter = _mm512_and_si512(w1, w2);
                lower = _mm512_extracti64x4_epi64(inter, 0);
                upper = _mm512_extracti64x4_epi64(inter, 1);
                CSA(twosB, ones, ones, lower, upper);
                CSA(foursB, twos, twos, twosA, twosB);
                CSA(eightsA, fours, fours, foursA, foursB);
                w1 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + k + 32));
                w2 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + k + 32));
                inter = _mm512_and_si512(w1, w2);
                lower = _mm512_extracti64x4_epi64(inter, 0);
                upper = _mm512_extracti64x4_epi64(inter, 1);
                CSA(twosA, ones, ones, lower, upper);
                w1 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + k + 40));
                w2 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + k + 40));
                inter = _mm512_and_si512(w1, w2);
                lower = _mm512_extracti64x4_epi64(inter, 0);
                upper = _mm512_extracti64x4_epi64(inter, 1);
                CSA(twosB, ones, ones, lower, upper);
                CSA(foursA, twos, twos, twosA, twosB);
                w1 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + k + 48));
                w2 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + k + 48));
                inter = _mm512_and_si512(w1, w2);
                lower = _mm512_extracti64x4_epi64(inter, 0);
                upper = _mm512_extracti64x4_epi64(inter, 1);
                CSA(twosA, ones, ones, lower, upper);
                w1 = _mm512_load_si512(
                    (__m512i *)(gt_tbl1 + i * words + k + 56));
                w2 = _mm512_load_si512(
                    (__m512i *)(gt_tbl2 + j * words + k + 56));
                inter = _mm512_and_si512(w1, w2);
                lower = _mm512_extracti64x4_epi64(inter, 0);
                upper = _mm512_extracti64x4_epi64(inter, 1);
                CSA(twosB, ones, ones, lower, upper);
                CSA(foursB, twos, twos, twosA, twosB);
                CSA(eightsB, fours, fours, foursA, foursB);
                CSA(sixteens, eights, eights, eightsA, eightsB);

                total = _mm256_add_epi64(total, popcount(sixteens));
            }

            total = _mm256_slli_epi64(total, 4); // * 16
            total = _mm256_add_epi64(
                total, _mm256_slli_epi64(popcount(eights), 3)); // += 8 * ...
            total = _mm256_add_epi64(
                total, _mm256_slli_epi64(popcount(fours), 2)); // += 4 * ...
            total = _mm256_add_epi64(
                total, _mm256_slli_epi64(popcount(twos), 1)); // += 2 * ...
            total = _mm256_add_epi64(total, popcount(ones));

            for (; k < words; k += 8) {
                __m512i w1 =
                    _mm512_load_si512((__m512i *)(gt_tbl1 + i * words + k));
                __m512i w2 =
                    _mm512_load_si512((__m512i *)(gt_tbl2 + j * words + k));
                __m512i inter = _mm512_and_si512(w1, w2);
                __m256i lower = _mm512_extracti64x4_epi64(inter, 0);
                __m256i upper = _mm512_extracti64x4_epi64(inter, 1);
                __m256i sum =
                    _mm256_add_epi64(popcount(upper), popcount(lower));
                total = _mm256_add_epi64(total, sum);
            }

            ct_tbl[i * 3 + j] =
                static_cast<uint64_t>(_mm256_extract_epi64(total, 0)) +
                static_cast<uint64_t>(_mm256_extract_epi64(total, 1)) +
                static_cast<uint64_t>(_mm256_extract_epi64(total, 2)) +
                static_cast<uint64_t>(_mm256_extract_epi64(total, 3));
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
    avx512bw_popcnt_256_harley_seal(t1.cases, t1.size, t2.cases, t1.cases_words,
                                    out.cases, out.size);
    avx512bw_popcnt_256_harley_seal(t1.ctrls, t1.size, t2.ctrls, t2.ctrls_words,
                                    out.ctrls, out.size);
}
