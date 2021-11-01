#include <bitset>
#include <fiuncho/Dataset.h>
#include <gtest/gtest.h>
#include <string>

std::string tped, tfam, rawfile, gametes_rawfile;

void checks(const Dataset<uint64_t> &d)
{
    EXPECT_EQ(10, d.snps);
    EXPECT_EQ(9000, d.cases);
    EXPECT_EQ(9800, d.ctrls);

#if ALIGN == 64
    EXPECT_EQ(144, d[0].cases_words);
    EXPECT_EQ(160, d[0].ctrls_words);
#elif ALIGN == 32
    EXPECT_EQ(144, d[0].cases_words);
    EXPECT_EQ(156, d[0].ctrls_words);
#else
    EXPECT_EQ(141, d[0].cases_words);
    EXPECT_EQ(154, d[0].ctrls_words);
#endif

    size_t words = d[0].cases_words;
    for (size_t i = 0; i < d.snps; i++) {
        size_t count = 0;
        for (size_t j = 0; j < 3; j++) {
            for (size_t w = 0; w < words; w++) {
                count += std::bitset<64>(d[i].cases[j * words + w]).count();
            }
        }
        EXPECT_EQ(d.cases, count);
    }
    words = d[0].ctrls_words;
    for (size_t i = 0; i < d.snps; i++) {
        size_t count = 0;
        for (size_t j = 0; j < 3; j++) {
            for (size_t w = 0; w < words; w++) {
                count += std::bitset<64>(d[i].ctrls[j * words + w]).count();
            }
        }
        EXPECT_EQ(d.ctrls, count);
    }
};

void equal(const Dataset<uint64_t> &d1, const Dataset<uint64_t> &d2)
{
    EXPECT_EQ(d1.snps, d2.snps);
    EXPECT_EQ(d1.cases, d2.cases);
    EXPECT_EQ(d1.ctrls, d2.ctrls);

    EXPECT_EQ(d1[0].cases_words, d2[0].cases_words);
    EXPECT_EQ(d1[0].ctrls_words, d2[0].ctrls_words);

    size_t words = d1[0].cases_words;
    for (size_t i = 0; i < d1.snps; i++) {
        for (size_t j = 0; j < 3; j++) {
            for (size_t w = 0; w < words; w++) {
                EXPECT_TRUE(d1[i].cases[j * words + w] ==
                            d2[i].cases[j * words + w]);
            }
        }
    }
    words = d1[0].ctrls_words;
    for (size_t i = 0; i < d1.snps; i++) {
        size_t count = 0;
        for (size_t j = 0; j < 3; j++) {
            for (size_t w = 0; w < words; w++) {
                EXPECT_EQ(d1[i].ctrls[j * words + w],
                          d2[i].ctrls[j * words + w]);
            }
        }
    }
}

namespace
{
TEST(Dataset, TPED)
{
    const Dataset<uint64_t> dataset = Dataset<uint64_t>::read({tped, tfam});
    checks(dataset);
};

TEST(Dataset, RAW)
{
    const Dataset<uint64_t> d_tped = Dataset<uint64_t>::read({tped, tfam});
    const Dataset<uint64_t> d_raw = Dataset<uint64_t>::read({rawfile});
    equal(d_tped, d_raw);
};

TEST(Dataset, GAMETES_RAW)
{
    const Dataset<uint64_t> d_tped = Dataset<uint64_t>::read({tped, tfam});
    const Dataset<uint64_t> d_gametes_raw =
        Dataset<uint64_t>::read({gametes_rawfile});
    equal(d_tped, d_gametes_raw);
};
} // namespace

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    assert(argc == 5); // gtest leaved unparsed arguments for you
    tped = argv[1];
    tfam = argv[2];
    rawfile = argv[3];
    gametes_rawfile = argv[4];
    return RUN_ALL_TESTS();
}