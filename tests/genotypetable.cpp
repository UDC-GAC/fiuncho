#include <bitset>
#include <fiuncho/ContingencyTable.h>
#include <fiuncho/GenotypeTable.h>
#include <fiuncho/Dataset.h>
#include <gtest/gtest.h>
#include <string>

std::string tped, tfam;

namespace
{
TEST(GenotypeTableTest, fill)
{
    const auto d = Dataset<uint64_t>::read({tped, tfam});
    const auto &t1 = d[0];
    const auto &t2 = d[1];
    GenotypeTable<uint64_t> result(2, t1.cases_words, t1.ctrls_words);

    EXPECT_EQ(9, result.size);
    EXPECT_EQ(t1.cases_words, result.cases_words);
    EXPECT_EQ(t1.ctrls_words, result.ctrls_words);

    GenotypeTable<uint64_t>::combine(t1, t2, result);

    for (size_t i = 0; i < result.size; i++) {
        for (size_t j = 0; j < result.cases_words; j++) {
            EXPECT_TRUE(result.cases[i * result.cases_words + j] ==
                        (t1.cases[(i / 3) * t1.cases_words + j] &
                         t2.cases[(i % 3) * t2.cases_words + j]));
        }
        for (size_t j = 0; j < result.ctrls_words; j++) {
            EXPECT_TRUE(result.ctrls[i * result.ctrls_words + j] ==
                        (t1.ctrls[(i / 3) * t1.ctrls_words + j] &
                         t2.ctrls[(i % 3) * t2.ctrls_words + j]));
        }
    }
}

TEST(GenotypeTableTest, popcnt)
{
    const auto d = Dataset<uint64_t>::read({tped, tfam});

    GenotypeTable<uint64_t> gt(2, d[0].cases_words, d[0].ctrls_words);
    GenotypeTable<uint64_t>::combine(d[0], d[1], gt);
    ContingencyTable<uint32_t> ct(2);
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
