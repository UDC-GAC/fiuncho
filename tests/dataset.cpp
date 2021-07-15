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
#include <fiuncho/dataset/Dataset.h>
#include <gtest/gtest.h>
#include <string>

std::string tped, tfam;

namespace
{
TEST(DatasetTest, Dataset)
{
#ifdef ALIGN
    const Dataset<uint64_t> dataset =
        Dataset<uint64_t>::read<ALIGN>(tped, tfam);
#else
    const Dataset<uint64_t> dataset = Dataset<uint64_t>::read(tped, tfam);
#endif

    EXPECT_EQ(10, dataset.snps);
    EXPECT_EQ(600, dataset.cases);
    EXPECT_EQ(1300, dataset.ctrls);
    EXPECT_EQ(1900, dataset.cases + dataset.ctrls);

#if ALIGN == 64
    EXPECT_EQ(16, dataset[0].cases_words);
    EXPECT_EQ(24, dataset[0].ctrls_words);
#elif ALIGN == 32
    EXPECT_EQ(12, dataset[0].cases_words);
    EXPECT_EQ(24, dataset[0].ctrls_words);
#else
    EXPECT_EQ(10, dataset[0].cases_words);
    EXPECT_EQ(21, dataset[0].ctrls_words);
#endif

    for (size_t i = 0; i < dataset.snps; i++) {
        size_t count = 0;
        for (size_t j = 0; j < 3; j++) {
            for (size_t w = 0; w < dataset[i].cases_words; w++) {
                count += std::bitset<64>(
                             dataset[i].cases[j * dataset[i].cases_words + w])
                             .count();
            }
        }
        EXPECT_EQ(dataset.cases, count);
    }

    for (size_t i = 0; i < dataset.snps; i++) {
        size_t count = 0;
        for (size_t j = 0; j < 3; j++) {
            for (size_t w = 0; w < dataset[i].ctrls_words; w++) {
                count += std::bitset<64>(
                             dataset[i].ctrls[j * dataset[i].ctrls_words + w])
                             .count();
            }
        }
        EXPECT_EQ(dataset.ctrls, count);
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