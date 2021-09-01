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

#include <bitset>
#include <cmath>
#include <fiuncho/algorithms/MutualInformation.h>

template <>
MutualInformation<float>::MutualInformation(unsigned int num_cases,
                                            unsigned int num_ctrls)
{
    inv_inds = 1.0 / (num_cases + num_ctrls);

    float p = num_cases * inv_inds;
    h_y = (-1.0) * p * logf(p);

    p = num_ctrls * inv_inds;
    h_y -= p * logf(p);
}

template <>
template <>
float MutualInformation<float>::compute<uint32_t>(
    const ContingencyTable<uint32_t> &table) const noexcept
{
    float p_case, p_ctrl, p_any, log_p_case, log_p_ctrl, log_p_any, h_x = 0.0,
                                                                    h_all = 0.0;
    // Constant references to objects are not enough for the compiler to
    // determine the number of loop iterations for vectorization, therefore they
    // have to be copied to a local constant. GCC also needs the "-ffast-math"
    // optimization parameter to vectorize this loop.
    const size_t iters = table.size;
    for (size_t i = 0; i < iters; i++) {
        p_case = table.cases[i] * inv_inds;
        h_all -= p_case * logf(p_case);
        p_ctrl = table.ctrls[i] * inv_inds;
        h_all -= p_ctrl * logf(p_ctrl);
        p_any = p_case + p_ctrl;
        h_x -= p_any * logf(p_any);
    }

    return h_x + h_y - h_all;
}