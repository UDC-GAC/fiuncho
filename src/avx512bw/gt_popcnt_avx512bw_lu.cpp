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

const __m512i lookup = _mm512_setr_epi64(
    0x0302020102010100llu, 0x0403030203020201llu, 0x0302020102010100llu,
    0x0403030203020201llu, 0x0302020102010100llu, 0x0403030203020201llu,
    0x0302020102010100llu, 0x0403030203020201llu);
const __m512i low_mask = _mm512_set1_epi8(0x0f);

/**
 * @brief Implements the combination of two genotype subtables and subsequent
 * contingency table computation. The function is a fork of the function
 * ``popcnt_AVX512BW_lookup_original`` from
 * https://github.com/WojciechMula/sse-popcount/tree/6feb3dba32c526b17de01e931c116900e3a23104,
 * modified to interleave the genotype table combination operations with the
 * `popcount` operations.
 */

inline void avx512bw_popcnt_lookup(const uint64_t *gt_tbl1,
                                   const size_t gt_size,
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
    avx512bw_popcnt_lookup(t1.cases, t1.size, t2.cases, t1.cases_words,
                           out.cases, out.size);
    avx512bw_popcnt_lookup(t1.ctrls, t1.size, t2.ctrls, t2.ctrls_words,
                           out.ctrls, out.size);
}
