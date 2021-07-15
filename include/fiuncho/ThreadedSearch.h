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
 * @file ThreadedSearch.h
 * @author Christian Ponte
 */

#ifndef FIUNCHO_THREADEDSEARCH_H
#define FIUNCHO_THREADEDSEARCH_H

#include <cmath>
#include <fiuncho/ContingencyTable.h>
#include <fiuncho/GenotypeTable.h>
#include <fiuncho/Search.h>
#include <fiuncho/algorithms/MutualInformation.h>
#include <fiuncho/dataset/Dataset.h>
#include <fiuncho/utils/MaxArray.h>
#include <iostream>
#include <pthread.h>
#include <thread>
#include <vector>

#define BLOCK_SIZE (int)(16384 / powf(3, args.order - 2))

/**
 * Epistasis search class that uses CPU multi-threading to complete the
 * search.
 */

class ThreadedSearch : public Search
{

    const unsigned int nthreads;

    class Args
    {
      public:
        const Dataset<uint64_t> &dataset;
        const unsigned short order;
        const Distribution<int> distribution;
        MaxArray<Result<int, float>> maxarray;
#ifdef BENCHMARK
        double elapsed_time;
        size_t combinations;
#endif

        Args(const Dataset<uint64_t> &dataset, const unsigned short order,
             const Distribution<int> &distribution, const size_t outputs)
            : dataset(dataset), order(order), distribution(distribution),
              maxarray(outputs)
        {
#ifdef BENCHMARK
            elapsed_time = 0;
            combinations = 0;
#endif
        }
    };

    static void thread_main(Args &args)
    {
#ifdef BENCHMARK
        struct timespec ts;
        if (clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts) == -1) {
            throw std::runtime_error("Error while CLOCK_THREAD_CPUTIME_ID");
        }
        double start_time = ts.tv_sec + ts.tv_nsec * 1e-9;
#endif
        if (args.order == 2) {
            search_order_2(args);
        } else {
            search_order_gt_2(args);
        }
#ifdef BENCHMARK
        if (clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts) == -1) {
            throw std::runtime_error("Error while CLOCK_THREAD_CPUTIME_ID");
        }
        args.elapsed_time = ts.tv_sec + ts.tv_nsec * 1e-9 - start_time;
#endif
    }

    static void search_order_2(Args &args)
    {
        int i, j, k;
        // Create the ContingencyTable vector, Result vector, and MI objects
        std::vector<ContingencyTable<uint32_t>> cts;
        std::vector<Result<int, float>> r(BLOCK_SIZE);
        cts.reserve(BLOCK_SIZE);
        for (i = 0; i < BLOCK_SIZE; ++i) {
            cts.emplace_back(2, args.dataset[0].cases_words,
                             args.dataset[0].ctrls_words);
            r[i].combination.resize(2);
        }
        MutualInformation<float> mi(args.dataset.cases, args.dataset.ctrls);
        // For each combination assigned by the distribution
        j = 0;
        for (auto c = args.distribution.begin(); c < args.distribution.end();
             ++c) {
            // Iterate over subsequent combinations
            for (i = c->back() + 1; i < (int)args.dataset.snps; ++i) {
                // If the block is full, compute all MI's
                if (j == BLOCK_SIZE) {
                    for (k = 0; k < BLOCK_SIZE; ++k) {
                        // Compute mutual information
                        r[k].val = mi.compute(cts[k]);
                        args.maxarray.add(r[k]);
                    }
#ifdef BENCHMARK
                    args.combinations += j;
#endif
                    j = 0;
                }
                r[j].combination[0] = c[0];
                r[j].combination[1] = i;
                // Fill contingency table
                GenotypeTable<uint64_t>::combine_and_popcnt(
                    args.dataset[c[0]], args.dataset[i], cts[j++]);
            }
        }
        // For each contingency table remaining in the block
        for (k = 0; k < j; ++k) {
            // Compute mutual information
            r[k].val = mi.compute(cts[k]);
            args.maxarray.add(r[k]);
        }
#ifdef BENCHMARK
        args.combinations += j;
#endif
    }

    static void search_order_gt_2(Args &args)
    {
        int i, j, k;
        // Allocate genotype tables of size < target interaction order
        std::vector<GenotypeTable<uint64_t>> gts;
        gts.reserve(args.order - 2);
        for (auto o = 2; o < args.order; ++o) {
            gts.emplace_back(o, args.dataset[0].cases_words,
                             args.dataset[0].ctrls_words);
        }
        // Create ContingencyTable vector, Result vector, and MI objects
        std::vector<ContingencyTable<uint32_t>> cts;
        std::vector<Result<int, float>> r(BLOCK_SIZE);
        cts.reserve(BLOCK_SIZE);
        for (i = 0; i < BLOCK_SIZE; ++i) {
            cts.emplace_back(args.order, args.dataset[0].cases_words,
                             args.dataset[0].ctrls_words);
            r[i].combination.resize(args.order);
        }
        MutualInformation<float> mi(args.dataset.cases, args.dataset.ctrls);
        // For each combination assigned by the distribution
        j = 0;
        for (auto c = args.distribution.begin(); c < args.distribution.end();
             ++c) {
            // Fill genotype tables
            GenotypeTable<uint64_t>::combine(args.dataset[c[0]],
                                             args.dataset[c[1]], gts[0]);
            for (auto i = 1; i < args.order - 2; ++i) {
                GenotypeTable<uint64_t>::combine(
                    gts[i - 1], args.dataset[c[i + 1]], gts[i]);
            }
            // Iterate over subsequent combinations
            for (i = c->back() + 1; i < (int)args.dataset.snps; ++i) {
                // If the block is full, compute all MI's
                if (j == BLOCK_SIZE) {
                    for (k = 0; k < BLOCK_SIZE; ++k) {
                        // Compute mutual information
                        r[k].val = mi.compute(cts[k]);
                        args.maxarray.add(r[k]);
                    }
#ifdef BENCHMARK
                    args.combinations += j;
#endif
                    j = 0;
                }
                memcpy(r[j].combination.data(), c->data(),
                       c->size() * sizeof(int));
                r[j].combination.back() = i;
                // Fill contingency table
                GenotypeTable<uint64_t>::combine_and_popcnt(
                    gts.back(), args.dataset[i], cts[j++]);
            }
        }
        // For each contingency table remaining in the block
        for (k = 0; k < j; ++k) {
            // Compute mutual information
            r[k].val = mi.compute(cts[k]);
            args.maxarray.add(r[k]);
        }
#ifdef BENCHMARK
        args.combinations += j;
#endif
    }

  public:
    /**
     * @name Constructors
     */
    //@{

    /**
     * Create a ThreadedSearch object.
     *
     * @param threads Number of threads to use during the search
     */

    ThreadedSearch(unsigned int threads) : nthreads(threads) {}

    //@}

    ~ThreadedSearch() {}

    /**
     * @name Methods
     */
    //@{

    std::vector<Result<int, float>> run(const Dataset<uint64_t> &dataset,
                                        const unsigned short order,
                                        const Distribution<int> &distribution,
                                        const unsigned int outputs)
    {
        // Spawn threads
        std::vector<Args> thread_args;
        std::vector<std::thread> threads;
        // Pre-reserve space to avoid reallocating the underlying array, which
        // results in an error since previous addresses are rendered incorrect
        thread_args.reserve(nthreads);
        threads.reserve(nthreads);
        for (unsigned int i = 0; i < nthreads; i++) {
            thread_args.emplace_back(dataset, order,
                                     distribution.layer(nthreads, i), outputs);
            threads.emplace_back(thread_main, std::ref(thread_args.back()));
        }

        std::vector<Result<int, float>> results;
        results.reserve(nthreads * outputs);
        // Wait for the completion of all threads
        for (unsigned int i = 0; i < threads.size(); i++) {
            threads[i].join();
            results.insert(
                results.end(), &thread_args[i].maxarray[0],
                &thread_args[i].maxarray[thread_args[i].maxarray.size()]);
#ifdef BENCHMARK
            // Print information
            std::cout << "Thread " << i << ": " << thread_args[i].elapsed_time
                      << "s, " << thread_args[i].combinations
                      << " combinations\n";
#endif
        }
        // Sort the auxiliar array and resize the result before returning
        std::sort(results.rbegin(), results.rend());
        if (results.size() > outputs) {
            results.resize(outputs);
        }
        return results;
    }

    //@}
};

#endif
