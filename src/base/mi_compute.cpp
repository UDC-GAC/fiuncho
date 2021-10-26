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
template <>
float MutualInformation<float>::compute<uint32_t>(
    const ContingencyTable<uint32_t> &table) const noexcept
{
    float p_case, p_ctrl, p_any, h_x = 0.0f, h_all = 0.0f;
    // Constant references to objects are not enough for the compiler to
    // determine the number of loop iterations for vectorization, therefore they
    // have to be copied to a local constant. GCC also needs the "-ffast-math"
    // optimization parameter to vectorize this loop.
    const size_t size = table.size;
    for (size_t i = 0; i < size; i++) {
        p_case = table.cases[i] * inv_inds;
        if (p_case != 0.0f) {
            h_all -= p_case * logf(p_case);
        }
        p_ctrl = table.ctrls[i] * inv_inds;
        if (p_ctrl != 0.0f) {
            h_all -= p_ctrl * logf(p_ctrl);
        }
        p_any = p_case + p_ctrl;
        if (p_any != 0.0f) {
            h_x -= p_any * logf(p_any);
        }
    }

    return h_x + h_y - h_all;
}
