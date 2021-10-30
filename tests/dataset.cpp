#include <bitset>
#include <fiuncho/Dataset.h>
#include <gtest/gtest.h>
#include <string>

std::string tped, tfam, rawfile;

auto checks = [](const Dataset<uint64_t> &dataset) {
    EXPECT_EQ(10, dataset.snps);
    EXPECT_EQ(9000, dataset.cases);
    EXPECT_EQ(9800, dataset.ctrls);

#if ALIGN == 64
    EXPECT_EQ(144, dataset[0].cases_words);
    EXPECT_EQ(160, dataset[0].ctrls_words);
#elif ALIGN == 32
    EXPECT_EQ(144, dataset[0].cases_words);
    EXPECT_EQ(156, dataset[0].ctrls_words);
#else
    EXPECT_EQ(141, dataset[0].cases_words);
    EXPECT_EQ(154, dataset[0].ctrls_words);
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
};

namespace
{
TEST(Dataset, TPED)
{
    const Dataset<uint64_t> dataset = Dataset<uint64_t>::read({tped, tfam});
    checks(dataset);
};

TEST(Dataset, RAW)
{
    const Dataset<uint64_t> dataset = Dataset<uint64_t>::read({rawfile});
    checks(dataset);
};
} // namespace

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    assert(argc == 4); // gtest leaved unparsed arguments for you
    tped = argv[1];
    tfam = argv[2];
    rawfile = argv[3];
    return RUN_ALL_TESTS();
}