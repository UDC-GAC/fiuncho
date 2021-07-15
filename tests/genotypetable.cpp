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

#include <bitset>
#include <fiuncho/ContingencyTable.h>
#include <fiuncho/GenotypeTable.h>
#include <gtest/gtest.h>

#ifdef ALIGN
alignas(ALIGN)
#endif
    uint64_t cases1[3][8] = {
        {0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0xffffffffffffffff, 0x0000000000000000, 0x0000000000000000,
         0x0000000000000000, 0x0000000000000000},
        {0xAAAAAAAAAAAAAAAA, 0xAAAAAAAAAAAAAAAA, 0x5555555555555555,
         0x5555555555555555, 0xAAAAAAAAAAAAAAAA, 0xAAAAAAAAAAAAAAAA,
         0x5555555555555555, 0x5555555555555555},
        {0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0xffffffffffffffff, 0xffffffffffffffff}};
#ifdef ALIGN
alignas(ALIGN)
#endif
    uint64_t cases2[3][8] = {
        {0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0xffffffffffffffff, 0x0000000000000000, 0x0000000000000000,
         0x0000000000000000, 0x0000000000000000},
        {0xAAAAAAAAAAAAAAAA, 0xAAAAAAAAAAAAAAAA, 0x5555555555555555,
         0x5555555555555555, 0xAAAAAAAAAAAAAAAA, 0xAAAAAAAAAAAAAAAA,
         0x5555555555555555, 0x5555555555555555},
        {0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0xffffffffffffffff, 0xffffffffffffffff}};
#ifdef ALIGN
alignas(ALIGN)
#endif
    uint64_t ctrls1[3][16] = {
        {0xAAAAAAAAAAAAAAAA, 0xAAAAAAAAAAAAAAAA, 0x5555555555555555,
         0x5555555555555555, 0xAAAAAAAAAAAAAAAA, 0xAAAAAAAAAAAAAAAA,
         0x5555555555555555, 0x5555555555555555, 0xAAAAAAAAAAAAAAAA,
         0xAAAAAAAAAAAAAAAA, 0x5555555555555555, 0x5555555555555555,
         0xAAAAAAAAAAAAAAAA, 0xAAAAAAAAAAAAAAAA, 0x5555555555555555,
         0x5555555555555555},
        {0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0xffffffffffffffff},
        {0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0xffffffffffffffff, 0x0000000000000000, 0x0000000000000000,
         0x0000000000000000, 0x0000000000000000, 0xffffffffffffffff,
         0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
         0x0000000000000000}};
#ifdef ALIGN
alignas(ALIGN)
#endif
    uint64_t ctrls2[3][16] = {
        {0xAAAAAAAAAAAAAAAA, 0xAAAAAAAAAAAAAAAA, 0x5555555555555555,
         0x5555555555555555, 0xAAAAAAAAAAAAAAAA, 0xAAAAAAAAAAAAAAAA,
         0x5555555555555555, 0x5555555555555555, 0xAAAAAAAAAAAAAAAA,
         0xAAAAAAAAAAAAAAAA, 0x5555555555555555, 0x5555555555555555,
         0xAAAAAAAAAAAAAAAA, 0xAAAAAAAAAAAAAAAA, 0x5555555555555555,
         0x5555555555555555},
        {0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0xffffffffffffffff},
        {0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0xffffffffffffffff, 0x0000000000000000, 0x0000000000000000,
         0x0000000000000000, 0x0000000000000000, 0xffffffffffffffff,
         0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
         0x0000000000000000}};

GenotypeTable<uint64_t> t1(cases1[0], 8, ctrls1[0], 16);
GenotypeTable<uint64_t> t2(cases2[0], 8, ctrls2[0], 16);

namespace
{
TEST(GenotypeTableTest, fill)
{
    uint64_t res_cases[9][8] = {
        {0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0xffffffffffffffff, 0x0000000000000000, 0x0000000000000000,
         0x0000000000000000, 0x0000000000000000},
        {0xAAAAAAAAAAAAAAAA, 0xAAAAAAAAAAAAAAAA, 0x5555555555555555,
         0x5555555555555555, 0x0000000000000000, 0x0000000000000000,
         0x0000000000000000, 0x0000000000000000},
        {0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0xffffffffffffffff, 0x0000000000000000, 0x0000000000000000,
         0x0000000000000000, 0x0000000000000000},
        {0xAAAAAAAAAAAAAAAA, 0xAAAAAAAAAAAAAAAA, 0x5555555555555555,
         0x5555555555555555, 0x0000000000000000, 0x0000000000000000,
         0x0000000000000000, 0x0000000000000000},
        {0xAAAAAAAAAAAAAAAA, 0xAAAAAAAAAAAAAAAA, 0x5555555555555555,
         0x5555555555555555, 0xAAAAAAAAAAAAAAAA, 0xAAAAAAAAAAAAAAAA,
         0x5555555555555555, 0x5555555555555555},
        {0xAAAAAAAAAAAAAAAA, 0xAAAAAAAAAAAAAAAA, 0x5555555555555555,
         0x5555555555555555, 0xAAAAAAAAAAAAAAAA, 0xAAAAAAAAAAAAAAAA,
         0x5555555555555555, 0x5555555555555555},
        {0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0xffffffffffffffff, 0x0000000000000000, 0x0000000000000000,
         0x0000000000000000, 0x0000000000000000},
        {0xAAAAAAAAAAAAAAAA, 0xAAAAAAAAAAAAAAAA, 0x5555555555555555,
         0x5555555555555555, 0xAAAAAAAAAAAAAAAA, 0xAAAAAAAAAAAAAAAA,
         0x5555555555555555, 0x5555555555555555},
        {0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0xffffffffffffffff, 0xffffffffffffffff}};
    uint64_t res_ctrls[9][16] = {
        {0xAAAAAAAAAAAAAAAA, 0xAAAAAAAAAAAAAAAA, 0x5555555555555555,
         0x5555555555555555, 0xAAAAAAAAAAAAAAAA, 0xAAAAAAAAAAAAAAAA,
         0x5555555555555555, 0x5555555555555555, 0xAAAAAAAAAAAAAAAA,
         0xAAAAAAAAAAAAAAAA, 0x5555555555555555, 0x5555555555555555,
         0xAAAAAAAAAAAAAAAA, 0xAAAAAAAAAAAAAAAA, 0x5555555555555555,
         0x5555555555555555},
        {0xAAAAAAAAAAAAAAAA, 0xAAAAAAAAAAAAAAAA, 0x5555555555555555,
         0x5555555555555555, 0xAAAAAAAAAAAAAAAA, 0xAAAAAAAAAAAAAAAA,
         0x5555555555555555, 0x5555555555555555, 0xAAAAAAAAAAAAAAAA,
         0xAAAAAAAAAAAAAAAA, 0x5555555555555555, 0x5555555555555555,
         0xAAAAAAAAAAAAAAAA, 0xAAAAAAAAAAAAAAAA, 0x5555555555555555,
         0x5555555555555555},
        {0xAAAAAAAAAAAAAAAA, 0xAAAAAAAAAAAAAAAA, 0x5555555555555555,
         0x5555555555555555, 0x0000000000000000, 0x0000000000000000,
         0x0000000000000000, 0x0000000000000000, 0xAAAAAAAAAAAAAAAA,
         0xAAAAAAAAAAAAAAAA, 0x5555555555555555, 0x5555555555555555,
         0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
         0x0000000000000000},
        {0xAAAAAAAAAAAAAAAA, 0xAAAAAAAAAAAAAAAA, 0x5555555555555555,
         0x5555555555555555, 0xAAAAAAAAAAAAAAAA, 0xAAAAAAAAAAAAAAAA,
         0x5555555555555555, 0x5555555555555555, 0xAAAAAAAAAAAAAAAA,
         0xAAAAAAAAAAAAAAAA, 0x5555555555555555, 0x5555555555555555,
         0xAAAAAAAAAAAAAAAA, 0xAAAAAAAAAAAAAAAA, 0x5555555555555555,
         0x5555555555555555},
        {0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0xffffffffffffffff},
        {0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0xffffffffffffffff, 0x0000000000000000, 0x0000000000000000,
         0x0000000000000000, 0x0000000000000000, 0xffffffffffffffff,
         0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
         0x0000000000000000},
        {0xAAAAAAAAAAAAAAAA, 0xAAAAAAAAAAAAAAAA, 0x5555555555555555,
         0x5555555555555555, 0x0000000000000000, 0x0000000000000000,
         0x0000000000000000, 0x0000000000000000, 0xAAAAAAAAAAAAAAAA,
         0xAAAAAAAAAAAAAAAA, 0x5555555555555555, 0x5555555555555555,
         0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
         0x0000000000000000},
        {0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0xffffffffffffffff, 0x0000000000000000, 0x0000000000000000,
         0x0000000000000000, 0x0000000000000000, 0xffffffffffffffff,
         0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
         0x0000000000000000},
        {0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0xffffffffffffffff, 0x0000000000000000, 0x0000000000000000,
         0x0000000000000000, 0x0000000000000000, 0xffffffffffffffff,
         0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
         0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
         0x0000000000000000}};

    GenotypeTable<uint64_t> result(2, 8, 16);

    EXPECT_EQ(9, result.size);
    EXPECT_EQ(8, result.cases_words);
    EXPECT_EQ(16, result.ctrls_words);

    GenotypeTable<uint64_t>::combine(t1, t2, result);

    for (size_t i = 0; i < result.size; i++) {
        for (size_t j = 0; j < result.cases_words; j++) {
            EXPECT_TRUE(result.cases[i * result.cases_words + j] ==
                        res_cases[i][j]);
        }
        for (size_t j = 0; j < result.ctrls_words; j++) {
            EXPECT_TRUE(result.ctrls[i * result.ctrls_words + j] ==
                        res_ctrls[i][j]);
        }
    }
}

TEST(GenotypeTableTest, popcnt)
{
    uint32_t popcnt_cases[9] = {256, 128, 256, 128, 256, 256, 256, 256, 512};
    uint64_t popcnt_ctrls[9] = {512, 512, 256, 512, 1024, 512, 256, 512, 512};

    ContingencyTable<uint32_t> ctable(2, 8, 16);

#ifdef ALIGN
    EXPECT_EQ(16, ctable.size);
#else
    EXPECT_EQ(9, ctable.size);
#endif
    EXPECT_EQ(8, ctable.cases_words);
    EXPECT_EQ(16, ctable.ctrls_words);

    GenotypeTable<uint64_t>::combine_and_popcnt(t1, t2, ctable);

    for (auto i = 0; i < 9; i++) {
        EXPECT_EQ(ctable.cases[i], popcnt_cases[i]);
    }
#ifdef ALIGN
    for (auto i = 9; i < 16; i++) {
        EXPECT_EQ(ctable.cases[i], 0);
    }
#endif

    for (auto i = 0; i < 9; i++) {
        EXPECT_EQ(ctable.ctrls[i], popcnt_ctrls[i]);
    }
#ifdef ALIGN
    for (auto i = 9; i < 16; i++) {
        EXPECT_EQ(ctable.ctrls[i], 0);
    }
#endif
}
} // namespace
