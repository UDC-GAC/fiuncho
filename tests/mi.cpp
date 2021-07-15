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

#include <gtest/gtest.h>
#include <fiuncho/ContingencyTable.h>
#include <fiuncho/algorithms/MutualInformation.h>

namespace
{
TEST(MI, compute)
{
    ContingencyTable<uint32_t> ctable(2, 0, 0);
#ifdef ALIGN
    for (size_t i = 9; i < ctable.size; i++) {
        ctable.cases[i] = 0;
        ctable.ctrls[i] = 0;
    }
#endif
    MutualInformation<float> mi(9, 18);

    for (size_t i = 0; i < 9; i++) {
        ctable.cases[i] = 1;
        ctable.ctrls[i] = 2;
    }
    EXPECT_NEAR(0, mi.compute(ctable), 1E-5);

    for (size_t i = 0; i < 9; i++) {
        ctable.cases[i] = 0;
        ctable.ctrls[i] = 0;
    }
    ctable.cases[0] = 9;
    ctable.ctrls[1] = 18;
    EXPECT_NEAR(0.6365141682948129, mi.compute(ctable), 1E-5);
}
} // namespace
