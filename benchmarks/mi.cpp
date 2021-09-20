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

/**
 * @file mi.cpp
 * @author Christian (christian.ponte@udc.es)
 * @brief Benchmarking program for the Mutual Information computation function
 * @date 20/09/2021
 */

#include "utils.h"
#include <fiuncho/ContingencyTable.h>
#include <fiuncho/GenotypeTable.h>
#include <fiuncho/algorithms/MutualInformation.h>
#include <fiuncho/dataset/Dataset.h>
#include <iostream>

void loop_mi(const int repetitions,
             std::vector<ContingencyTable<uint32_t>> &ctables,
             const MutualInformation<float> &mi)
{
    for (auto reps = 0; reps < repetitions; reps++) {
        for (auto &ctable : ctables) {
            mi.compute(ctable);
        }
    }
}

/**
 * @brief Main benchmark function. It precomputes a list of contingency tables
 * and measures the elapsed time during the calculation of their Mutual
 * Information for a given number of repetitions. The thread is pinned to the
 * indicated CPU core, and every thread is synchronized at the start.
 *
 * @param tped TPED file path
 * @param tfam TFAM file path
 * @param order Order of the input contingency tables to the MI function
 * @param affinity List of CPU cores to be used
 * @param repetitions Number of repetitions to loop during the measurement
 * @return A vector containing the elapsed times for each thread
 */

std::vector<double> bench_mi(const std::string tped, const std::string tfam,
                             const unsigned short order,
                             const std::vector<int> &affinity,
                             const int repetitions)
{
#ifdef ALIGN
    const auto dataset = Dataset<uint64_t>::read<ALIGN>(tped, tfam);
#else
    const auto dataset = Dataset<uint64_t>::read(tped, tfam);
#endif
    const size_t snp_count = dataset.snps, cases_words = dataset[0].cases_words,
                 ctrls_words = dataset[0].ctrls_words;
    // Initialize contingency tables outside of the main loop
    std::vector<GenotypeTable<uint64_t>> gtables;
    gtables.reserve(order - 2);
    std::vector<ContingencyTable<uint32_t>> ctables;
    ctables.reserve(snp_count - (order - 1));
    if (order > 2) {
        gtables.emplace_back(2, cases_words, ctrls_words);
        GenotypeTable<uint64_t>::combine(dataset[0], dataset[1], gtables[0]);
        for (auto o = 3; o < order; o++) {
            gtables.emplace_back(o, cases_words, ctrls_words);
            GenotypeTable<uint64_t>::combine(gtables[o - 3], dataset[o - 1],
                                             gtables[o - 2]);
        }
    }
    for (size_t snp = order - 1; snp < snp_count; snp++) {
        ctables.emplace_back(order, cases_words, ctrls_words);
        GenotypeTable<uint64_t>::combine_and_popcnt(
            (order > 2) ? gtables.back() : dataset[0], dataset[snp],
            ctables.back());
    }
    MutualInformation<float> mi(dataset.cases, dataset.ctrls);
    // Run benchmark function in parallel
    return multithread_pinned_time(affinity, loop_mi, repetitions, ctables, mi);
}

int main(int argc, char *argv[])
{
    if (argc != 6) {
        std::cout << argv[0] << " <THREADS> <ORDER> <REPETITIONS> <TPED> <TFAM>"
                  << std::endl;
        return 0;
    }
    // Arguments
    std::vector<int> affinity = split_into_ints(std::string(argv[1]), ',');
    const unsigned short order = atoi(argv[2]);
    const int repetitions = atoi(argv[3]);
    const std::string tped = argv[4], tfam = argv[5];
    // Run benchmark
    auto times = bench_mi(tped, tfam, order, affinity, repetitions);
    // Print times
    for (size_t i = 0; i < times.size() - 1; ++i) {
        std::cout << times[i] << ',';
    }
    std::cout << times.back() << '\n';

    return 0;
}
