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

/**
 * @brief Implements the combination of two genotype subtables and subsequent
 * contingency table computation. The function is a fork of the function
 * ``builtin_popcnt_movdq`` from
 * https://github.com/WojciechMula/sse-popcount/tree/6feb3dba32c526b17de01e931c116900e3a23104
 * modified to interleave the genotype table combination operations with the
 * `popcount` operations.
 */

inline void avx2_popcnt_64_native_movdq(const uint64_t *gt_tbl1,
                                        const size_t gt_size,
                                        const uint64_t *gt_tbl2,
                                        const size_t words, uint32_t *ct_tbl,
                                        const size_t ct_size)
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
                const __m128i lo = _mm256_extracti128_si256(inters, 0);
                uint64_t lower64 = _mm_cvtsi128_si64(lo);
                cnt += _mm_popcnt_u64(lower64);
                __m128i temp2 = (__m128i)_mm_movehl_ps((__m128)lo, (__m128)lo);
                uint64_t upper64 = _mm_cvtsi128_si64(temp2);
                cnt += _mm_popcnt_u64(upper64);
                const __m128i hi = _mm256_extracti128_si256(inters, 1);
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
    avx2_popcnt_64_native_movdq(t1.cases, t1.size, t2.cases, t1.cases_words,
                                out.cases, out.size);
    avx2_popcnt_64_native_movdq(t1.ctrls, t1.size, t2.ctrls, t2.ctrls_words,
                                out.ctrls, out.size);
}
