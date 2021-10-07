/*
 * This file is part of Fiuncho.
 *
 * Fiuncho is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Fiuncho is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Fiuncho. If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * @file cpu/MutualInformation.cpp
 * @author Jorge González
 * @author Christian Ponte
 * @date 1 October 2018
 *
 * @brief MutualInformation class members implementation.
 */

#include <cmath>
#include <fiuncho/algorithms/MutualInformation.h>
#include <immintrin.h>

#if !defined(__INTEL_COMPILER) && !defined(__INTEL_LLVM_COMPILER)
extern "C" {
__m256 _ZGVdN8v_logf(__m256 x);
}
inline __m256 _mm256_log_ps(__m256 x) noexcept { return _ZGVdN8v_logf(x); }
#endif

template <>
MutualInformation<float>::MutualInformation(unsigned int num_cases,
                                            unsigned int num_ctrls)
{
    inv_inds = 1.0 / (num_cases + num_ctrls);

    float p = num_cases * inv_inds;
    h_y = p * logf(p);

    p = num_ctrls * inv_inds;
    h_y += p * logf(p);
}

template <>
template <>
float MutualInformation<float>::compute<uint32_t>(
    const ContingencyTable<uint32_t> &table) const noexcept
{
    const __m256 ones = _mm256_set1_ps(1.0), ii = _mm256_set1_ps(inv_inds);

    __m256 y0, y1, y2, y3, y4, y5;

    __m256 h_x = _mm256_setzero_ps(), h_all = _mm256_setzero_ps();
    for (size_t i = 0; i < table.size; i += 8) {
        y0 = _mm256_cvtepi32_ps(
            _mm256_load_si256((const __m256i *)(table.cases + i)));
        y3 = _mm256_mul_ps(y0, ii);
        // Identify values different from 0
        y1 = _mm256_cmp_ps(y0, _mm256_setzero_ps(), _CMP_NEQ_OQ);
        // Replace 0's with 1's
        y4 = _mm256_log_ps(_mm256_blendv_ps(ones, y3, y1));
        h_all = _mm256_fmadd_ps(y3, y4, h_all);

        y0 = _mm256_cvtepi32_ps(
            _mm256_load_si256((const __m256i *)(table.ctrls + i)));
        y4 = _mm256_mul_ps(y0, ii);
        // Identify values different from 0
        y2 = _mm256_cmp_ps(y0, _mm256_setzero_ps(), _CMP_NEQ_OQ);
        // Replace 0's with 1's
        y5 = _mm256_log_ps(_mm256_blendv_ps(ones, y4, y2));
        h_all = _mm256_fmadd_ps(y4, y5, h_all);
        y5 = _mm256_add_ps(y3, y4);
        // Merge previous masks
        y1 = _mm256_or_ps(y1, y2);
        // Replace 0's with 1's
        y3 = _mm256_log_ps(_mm256_blendv_ps(ones, y5, y1));
        h_x = _mm256_fmadd_ps(y5, y3, h_x);
    }

    __m256 h_sum = _mm256_hadd_ps(h_all, h_x);
    return (h_sum[0] + h_sum[1] + h_sum[4] + h_sum[5]) - h_y -
           (h_sum[2] + h_sum[3] + h_sum[6] + h_sum[7]);
}
