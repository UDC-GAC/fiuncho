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
__m512 _ZGVeN16v_logf(__m512 x);
}
inline __m512 _mm512_log_ps(__m512 x) noexcept { return _ZGVeN16v_logf(x); }
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
    const __m512 ones = _mm512_set1_ps(1.0), ii = _mm512_set1_ps(inv_inds);
    __m512 h_x = _mm512_setzero_ps(), h_all = _mm512_setzero_ps();

    __m512i z0;
    __m512 z1, z2, z3, z4;
    __mmask16 mask1, mask2, mask3;
    for (size_t i = 0; i < table.size; i += 16) {
        z0 = _mm512_load_si512(table.cases + i);
        mask1 =
            _mm512_cmp_epi32_mask(z0, _mm512_setzero_si512(), _MM_CMPINT_NE);
        z1 = _mm512_cvtepi32_ps(z0);
        z2 = _mm512_mul_ps(z1, ii);
#ifdef __INTEL_COMPILER
        z3 = _mm512_mask_log_ps(_mm512_setzero_ps(), mask1, z2);
#else
        z3 = _mm512_log_ps(_mm512_mask_blend_ps(mask1, ones, z2));
#endif
        h_all = _mm512_fmadd_ps(z2, z3, h_all);

        z0 = _mm512_load_si512(table.ctrls + i);
        mask2 =
            _mm512_cmp_epi32_mask(z0, _mm512_setzero_si512(), _MM_CMPINT_NE);
        z1 = _mm512_cvtepi32_ps(z0);
        z3 = _mm512_mul_ps(z1, ii);
#ifdef __INTEL_COMPILER
        z4 = _mm512_mask_log_ps(_mm512_setzero_ps(), mask2, z3);
#else
        z4 = _mm512_log_ps(_mm512_mask_blend_ps(mask2, ones, z3));
#endif
        h_all = _mm512_fmadd_ps(z3, z4, h_all);

        mask3 = _kor_mask16(mask1, mask2);
        z1 = _mm512_add_ps(z2, z3);
#ifdef __INTEL_COMPILER
        z2 = _mm512_mask_log_ps(_mm512_setzero_ps(), mask3, z1);
#else
        z2 = _mm512_log_ps(_mm512_mask_blend_ps(mask3, ones, z1));
#endif
        h_x = _mm512_fmadd_ps(z1, z2, h_x);
    }
    return _mm512_reduce_add_ps(h_all) - _mm512_reduce_add_ps(h_x) - h_y;
}
