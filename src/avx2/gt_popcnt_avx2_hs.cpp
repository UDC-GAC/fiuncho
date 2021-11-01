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
 * https://github.com/WojciechMula/sse-popcount/tree/6feb3dba32c526b17de01e931c116900e3a23104
 * modified to interleave the genotype table combination operations with the
 * `popcount` operations.
 */

inline void avx2_popcnt_harley_seal(const uint64_t *gt_tbl1,
                                    const size_t gt_size,
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
                __m256i w11 =
                    _mm256_load_si256((__m256i *)(gt_tbl1 + i * words + k + 0));
                __m256i w21 =
                    _mm256_load_si256((__m256i *)(gt_tbl2 + j * words + k + 0));
                __m256i i1 = _mm256_and_si256(w11, w21);
                __m256i w12 =
                    _mm256_load_si256((__m256i *)(gt_tbl1 + i * words + k + 4));
                __m256i w22 =
                    _mm256_load_si256((__m256i *)(gt_tbl2 + j * words + k + 4));
                __m256i i2 = _mm256_and_si256(w12, w22);
                CSA(twosA, ones, ones, i1, i2);
                w11 =
                    _mm256_load_si256((__m256i *)(gt_tbl1 + i * words + k + 8));
                w21 =
                    _mm256_load_si256((__m256i *)(gt_tbl2 + j * words + k + 8));
                i1 = _mm256_and_si256(w11, w21);
                w12 = _mm256_load_si256(
                    (__m256i *)(gt_tbl1 + i * words + k + 12));
                w22 = _mm256_load_si256(
                    (__m256i *)(gt_tbl2 + j * words + k + 12));
                i2 = _mm256_and_si256(w12, w22);
                CSA(twosB, ones, ones, i1, i2);
                CSA(foursA, twos, twos, twosA, twosB);
                w11 = _mm256_load_si256(
                    (__m256i *)(gt_tbl1 + i * words + k + 16));
                w21 = _mm256_load_si256(
                    (__m256i *)(gt_tbl2 + j * words + k + 16));
                i1 = _mm256_and_si256(w11, w21);
                w12 = _mm256_load_si256(
                    (__m256i *)(gt_tbl1 + i * words + k + 20));
                w22 = _mm256_load_si256(
                    (__m256i *)(gt_tbl2 + j * words + k + 20));
                i2 = _mm256_and_si256(w12, w22);
                CSA(twosA, ones, ones, i1, i2);
                w11 = _mm256_load_si256(
                    (__m256i *)(gt_tbl1 + i * words + k + 24));
                w21 = _mm256_load_si256(
                    (__m256i *)(gt_tbl2 + j * words + k + 24));
                i1 = _mm256_and_si256(w11, w21);
                w12 = _mm256_load_si256(
                    (__m256i *)(gt_tbl1 + i * words + k + 28));
                w22 = _mm256_load_si256(
                    (__m256i *)(gt_tbl2 + j * words + k + 28));
                i2 = _mm256_and_si256(w12, w22);
                CSA(twosB, ones, ones, i1, i2);
                CSA(foursB, twos, twos, twosA, twosB);
                CSA(eightsA, fours, fours, foursA, foursB);
                w11 = _mm256_load_si256(
                    (__m256i *)(gt_tbl1 + i * words + k + 32));
                w21 = _mm256_load_si256(
                    (__m256i *)(gt_tbl2 + j * words + k + 32));
                i1 = _mm256_and_si256(w11, w21);
                w12 = _mm256_load_si256(
                    (__m256i *)(gt_tbl1 + i * words + k + 36));
                w22 = _mm256_load_si256(
                    (__m256i *)(gt_tbl2 + j * words + k + 36));
                i2 = _mm256_and_si256(w12, w22);
                CSA(twosA, ones, ones, i1, i2);
                w11 = _mm256_load_si256(
                    (__m256i *)(gt_tbl1 + i * words + k + 40));
                w21 = _mm256_load_si256(
                    (__m256i *)(gt_tbl2 + j * words + k + 40));
                i1 = _mm256_and_si256(w11, w21);
                w12 = _mm256_load_si256(
                    (__m256i *)(gt_tbl1 + i * words + k + 44));
                w22 = _mm256_load_si256(
                    (__m256i *)(gt_tbl2 + j * words + k + 44));
                i2 = _mm256_and_si256(w12, w22);
                CSA(twosB, ones, ones, i1, i2);
                CSA(foursA, twos, twos, twosA, twosB);
                w11 = _mm256_load_si256(
                    (__m256i *)(gt_tbl1 + i * words + k + 48));
                w21 = _mm256_load_si256(
                    (__m256i *)(gt_tbl2 + j * words + k + 48));
                i1 = _mm256_and_si256(w11, w21);
                w12 = _mm256_load_si256(
                    (__m256i *)(gt_tbl1 + i * words + k + 52));
                w22 = _mm256_load_si256(
                    (__m256i *)(gt_tbl2 + j * words + k + 52));
                i2 = _mm256_and_si256(w12, w22);
                CSA(twosA, ones, ones, i1, i2);
                w11 = _mm256_load_si256(
                    (__m256i *)(gt_tbl1 + i * words + k + 56));
                w21 = _mm256_load_si256(
                    (__m256i *)(gt_tbl2 + j * words + k + 56));
                i1 = _mm256_and_si256(w11, w21);
                w12 = _mm256_load_si256(
                    (__m256i *)(gt_tbl1 + i * words + k + 60));
                w22 = _mm256_load_si256(
                    (__m256i *)(gt_tbl2 + j * words + k + 60));
                i2 = _mm256_and_si256(w12, w22);
                CSA(twosB, ones, ones, i1, i2);
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

            for (; k < words; k += 4) {
                __m256i w1 =
                    _mm256_load_si256((__m256i *)(gt_tbl1 + i * words + k));
                __m256i w2 =
                    _mm256_load_si256((__m256i *)(gt_tbl2 + j * words + k));
                __m256i inters = _mm256_and_si256(w1, w2);
                total = _mm256_add_epi64(total, popcount(inters));
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
    avx2_popcnt_harley_seal(t1.cases, t1.size, t2.cases, t1.cases_words,
                            out.cases, out.size);
    avx2_popcnt_harley_seal(t1.ctrls, t1.size, t2.ctrls, t2.ctrls_words,
                            out.ctrls, out.size);
}
