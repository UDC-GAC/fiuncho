#include "utils.h"
#include <fiuncho/Distribution.h>
#include <fiuncho/Search.h>
#include <fiuncho/ThreadedSearch.h>
#include <fiuncho/Dataset.h>
#include <gtest/gtest.h>

std::string tped, tfam;

namespace
{
TEST(ThreadedSearchTest, Main)
{
    // Run ThreadedSearch
    const auto dataset = Dataset<uint64_t>::read({tped, tfam});

    std::vector<int> thread_count_vector{1, 32};
    for (auto t : thread_count_vector) {
        ThreadedSearch search(t);
        for (auto o = 2; o < 5; o++) {
            Distribution<int> distribution(dataset.snps, o - 1, 1, 0);
            auto result = search.run(dataset, o, distribution, 90);
            EXPECT_GT(result.size(), 0);
            if (o > 2) {
                EXPECT_EQ(result.size(), 90);
            }
            EXPECT_EQ(result[0].combination.size(), o);
            EXPECT_FALSE(has_repeated_elements(result));
            EXPECT_TRUE(ascending_combinations(result));
            if (o == 3) {
                EXPECT_TRUE(matches_mpi3snp_output(result));
            }
        }
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