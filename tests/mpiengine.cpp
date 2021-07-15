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

#include "utils.h"
#include <fiuncho/MPIEngine.h>
#include <fiuncho/ThreadedSearch.h>
#include <fiuncho/dataset/Dataset.h>
#include <gtest/gtest.h>

std::string tped, tfam;

namespace
{
TEST(MPIEngineTest, General)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPIEngine engine;
    for (auto o = 2; o < 5; o++) {
        auto results = engine.run<ThreadedSearch>(tped, tfam, o, 100, 4);
        if (rank == 0){
            EXPECT_GT(results.size(), 0);
            if (o > 2){
                EXPECT_EQ(results.size(), 100);
            }
            EXPECT_EQ(results[0].combination.size(), o);
            EXPECT_FALSE(has_repeated_elements(results));
            EXPECT_TRUE(ascending_combinations(results));
            if (o == 3) {
                EXPECT_TRUE(matches_mpi3snp_output(results));
            }
        }
    }
}
} // namespace

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleTest(&argc, argv);

    // Delete listeners for all processes minus the first
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ::testing::TestEventListeners& listeners =
    ::testing::UnitTest::GetInstance()->listeners();
    if (rank != 0) {
        delete listeners.Release(listeners.default_result_printer());
    }
    assert(argc == 3); // gtest leaved unparsed arguments for you
    tped = argv[1];
    tfam = argv[2];
    auto returncode = RUN_ALL_TESTS();
    MPI_Finalize();
    return returncode;
}