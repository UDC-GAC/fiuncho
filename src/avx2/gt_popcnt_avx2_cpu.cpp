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
#include <x86intrin.h>

inline void iter(const uint64_t *ptr1, const uint64_t *ptr2, uint32_t &result)
{
    __m256i w1 = _mm256_load_si256((__m256i *)ptr1);
    __m256i w2 = _mm256_load_si256((__m256i *)ptr2);
    __m256i inter = _mm256_and_si256(w1, w2);
    result += _popcnt64(_mm256_extract_epi64(inter, 0));
    result += _popcnt64(_mm256_extract_epi64(inter, 1));
    result += _popcnt64(_mm256_extract_epi64(inter, 2));
    result += _popcnt64(_mm256_extract_epi64(inter, 3));
}

/**
 * @brief Implements the combination of two genotype subtables and subsequent
 * contingency table computation. The function is a fork of the function
 * ``popcnt_AVX2_and_cpu`` from
 * https://github.com/WojciechMula/sse-popcount/tree/6feb3dba32c526b17de01e931c116900e3a23104,
 * modified to interleave the genotype table combination operations with the
 * `popcount` operations.
 */

inline void avx2_popcnt_cpu(const uint64_t *gt_tbl1, const size_t gt_size,
                            const uint64_t *gt_tbl2, const size_t words,
                            uint32_t *ct_tbl, const size_t ct_size)
{
    size_t i, j, k;
    for (i = 0; i < gt_size; ++i) {
        for (j = 0; j < 3; ++j) {
            uint32_t result = 0;
            for (k = 0; k + 16 <= words; k += 16) {
                iter(gt_tbl1 + i * words + k + 0, gt_tbl2 + j * words + k + 0,
                     result);
                iter(gt_tbl1 + i * words + k + 4, gt_tbl2 + j * words + k + 4,
                     result);
                iter(gt_tbl1 + i * words + k + 8, gt_tbl2 + j * words + k + 8,
                     result);
                iter(gt_tbl1 + i * words + k + 12, gt_tbl2 + j * words + k + 12,
                     result);
            }
            for (; k < words; k += 4) {
                iter(gt_tbl1 + i * words + k, gt_tbl2 + j * words + k, result);
            }
            ct_tbl[i * 3 + j] = result;
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
    avx2_popcnt_cpu(t1.cases, t1.size, t2.cases, t1.cases_words, out.cases,
                    out.size);
    avx2_popcnt_cpu(t1.ctrls, t1.size, t2.ctrls, t1.ctrls_words, out.ctrls,
                    out.size);
}
