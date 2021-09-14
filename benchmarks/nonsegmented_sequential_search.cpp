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
 * @file singlethreadsearch.cpp
 * @author Christian Ponte (christian.ponte@udc.es)
 * @brief Benchmark program for running single-thread epistasis searches
 * replicated on multiple (pinned) cores.
 * @date 2021-09-10
 */

#include "utils.h"
#include <cmath>
#include <cstdint>
#include <fiuncho/GenotypeTable.h>
#include <fiuncho/algorithms/MutualInformation.h>
#include <fiuncho/dataset/Dataset.h>
#include <fiuncho/utils/MaxArray.h>
#include <fiuncho/utils/Result.h>
#include <fiuncho/utils/StaticStack.h>
#include <iostream>

/**
 * @brief Single thread epistasis search. The method implements a depth-first
 * strategy using a static stack.
 *
 * @param dataset Input data
 * @param order Order of the search
 * @param[out] maxarray List of combinations with the largest MI value
 */

void stack_search(const Dataset<uint64_t> &dataset, const unsigned short order,
                  MaxArray<Result<int, float>> &maxarray)
{
    // Auxiliary definitions
    const size_t cases_words = dataset[0].cases_words,
                 ctrls_words = dataset[0].ctrls_words;
    const int snp_count = dataset.snps;
    // Define static structure for storing index combinations
    typedef struct {
        unsigned short size;
        int pos[];
    } Combination;
    const size_t item_size = sizeof(Combination) + order * sizeof(int);
    // Use a stack to explore combination space
    StaticStack<Combination> stack(item_size, snp_count * order - 2);
    Combination *cbuffer = (Combination *)new char[item_size];
    // Auxiliary bit tables for combinations sized under the target order
    std::vector<GenotypeTable<uint64_t>> btables;
    for (auto o = 2; o < order; o++) {
        btables.emplace_back(o, cases_words, ctrls_words);
    }
    // Vector of contingency tables (and their SNPs) for block processing
    ContingencyTable<uint32_t> ctable(order, cases_words, ctrls_words);
    MutualInformation<float> mi(dataset.cases, dataset.ctrls);
    Result<int, float> r;
    r.combination.resize(order);

    int i, j;
    for (i = 0; i < snp_count; ++i) {
        // Add all pairs starting with i to the stack
        cbuffer->size = 2;
        cbuffer->pos[0] = i;
        for (j = i + 1; j < snp_count; ++j) {
            cbuffer->pos[1] = j;
            stack.push(*cbuffer);
        }
        while (!stack.empty()) {
            // Process the combination from the top of the stack
            stack.pop(*cbuffer);
            const auto &last = cbuffer->pos[cbuffer->size - 1];
            // If the combination size is the same as the search order
            if (cbuffer->size == order) {
                // Compute the contingency table
                GenotypeTable<uint64_t>::combine_and_popcnt(
                    (order == 2) ? dataset[cbuffer->pos[0]] : btables.back(),
                    dataset[last], ctable);
                // Compute its MI value
                r.val = mi.compute(ctable);
                // Add to the array of results
                memcpy(r.combination.data(), cbuffer->pos, order * sizeof(int));
                maxarray.add(r);
            } else {
                // Compute the genotype table
                GenotypeTable<uint64_t>::combine(
                    (cbuffer->size == 2) ? dataset[cbuffer->pos[0]]
                                         : btables[cbuffer->size - 3],
                    dataset[last], btables[cbuffer->size - 2]);
                // Add subsequent combinations to top of the stack
                cbuffer->size += 1;
                for (j = last + 1; j < snp_count; ++j) {
                    cbuffer->pos[cbuffer->size - 1] = j;
                    stack.push(*cbuffer);
                }
            }
        }
    }
    free(cbuffer);
}

/**
 * @brief Meassure a single thread epistasis search. The thread is pinned to the
 * cpu core selected.
 *
 * @param tped Input TPED file
 * @param tfam Input TFAM file
 * @param order Order of the search
 * @param affinity Cpu core to which the thread will be pinned
 * @param[out] elapsed_time Elapsed time during the search
 */

void search(const std::string tped, const std::string tfam,
            const unsigned short order, const int affinity,
            double &elapsed_time)
{
    // Prepare input and output variables
#ifdef ALIGN
    const auto dataset = Dataset<uint64_t>::read<ALIGN>(tped, tfam);
#else
    const auto dataset = Dataset<uint64_t>::read(tped, tfam);
#endif
    MaxArray<Result<int, float>> result(10);
    // Measure search time
    elapsed_time = pinned_time(affinity, stack_search, dataset, order, result);
}

int main(int argc, char *argv[])
{
    if (argc != 5) {
        std::cout << argv[0] << " <THREADS> <ORDER> <TPED> <TFAM>" << std::endl;
        return 0;
    }

    // Initialization
    // Arguments
    std::vector<int> affinity = split_into_ints(std::string(argv[1]), ',');
    int thread_count = affinity.size();
    const unsigned short order = atoi(argv[2]);
    const std::string tped = argv[3], tfam = argv[4];
    // Variables
    pthread_barrier_init(&barrier, NULL, thread_count);
    std::vector<std::thread> threads;
    std::vector<double> times;
    times.resize(thread_count);

    // Spawn thread_count - 1 threads
    for (size_t i = 1; i < affinity.size(); i++) {
        threads.emplace_back(search, tped, tfam, order, affinity[i],
                             std::ref(times[i]));
    }
    // Also use current thread
    search(tped, tfam, order, affinity[0], times[0]);

    // Finalization
    // Wait for completion
    for (auto &thread : threads) {
        thread.join();
    }
    // Print times
    for (auto i = 0; i < thread_count - 1; i++) {
        std::cout << times[i] << ',';
    }
    std::cout << times[thread_count - 1] << '\n';

    return 0;
}
