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
#include <fiuncho/dataset/Dataset.h>
#include <gtest/gtest.h>
#include <string>

std::string tped, tfam;

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
#ifdef ALIGN
    const auto d = Dataset<uint64_t>::read<ALIGN>(tped, tfam);
#else
    const auto d = Dataset<uint64_t>::read(tped, tfam);
#endif

    GenotypeTable<uint64_t> gt(2, d[0].cases_words, d[0].ctrls_words);
    GenotypeTable<uint64_t>::combine(d[0], d[1], gt);
    ContingencyTable<uint32_t> ct(2, d[0].cases_words, d[0].ctrls_words);
    GenotypeTable<uint64_t>::combine_and_popcnt(d[0], d[1], ct);

    size_t i, j, k;
    const size_t cases_words = gt.cases_words, ctrls_words = gt.ctrls_words;
    for (i = 0; i < gt.size; ++i) {
        uint32_t cnt = 0;
        for (j = 0; j < cases_words; ++j) {
            cnt += std::bitset<64>(gt.cases[i * cases_words + j]).count();
        }
        EXPECT_EQ(cnt, ct.cases[i]);
        cnt = 0;
        for (j = 0; j < ctrls_words; ++j) {
            cnt += std::bitset<64>(gt.ctrls[i * ctrls_words + j]).count();
        }
        EXPECT_EQ(cnt, ct.ctrls[i]);
    }
    for (; i < ct.size; ++i) {
        EXPECT_EQ(0, ct.cases[i]);
        EXPECT_EQ(0, ct.ctrls[i]);
    }
}
} // namespace

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    assert(argc == 3); // gtest leaved unparsed arguments for you
    tped = argv[1];
    tfam = argv[2];
    return RUN_ALL_TESTS();
}
