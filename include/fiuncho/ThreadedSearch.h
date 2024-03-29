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
#include <fiuncho/Dataset.h>
#include <fiuncho/utils/MaxArray.h>
#include <iostream>
#include <pthread.h>
#include <thread>
#include <vector>

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

#ifdef SEGMENT_SIZE
        segmented_search(args, SEGMENT_SIZE);
#else
        nonsegmented_search(args);
#endif

#ifdef BENCHMARK
        if (clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts) == -1) {
            throw std::runtime_error("Error while CLOCK_THREAD_CPUTIME_ID");
        }
        args.elapsed_time = ts.tv_sec + ts.tv_nsec * 1e-9 - start_time;
#endif
    }

    static void segmented_search(Args &args, const size_t BASE_B)
    {
        const size_t B = BASE_B / (int)powf(3, args.order);
        int i, j, k;
        // Allocate genotype tables of size < target interaction order
        std::vector<GenotypeTable<uint64_t>> gts;
        gts.reserve(args.order - 2);
        for (auto o = 2; o < args.order; ++o) {
            gts.emplace_back(o, args.dataset[0].cases_words,
                             args.dataset[0].ctrls_words);
        }
        // Create ContingencyTable vector, Result vector, and MI objects
        auto cts_array = ContingencyTable<uint32_t>::make_array(B, args.order);
        auto cts = cts_array.get();
        std::vector<Result<int, float>> r(B);
        for (i = 0; i < B; ++i) {
            r[i].combination.resize(args.order);
        }
        MutualInformation<float> mi(args.dataset.cases, args.dataset.ctrls);
        // For each combination assigned by the distribution
        j = 0;
        for (auto c = args.distribution.begin(); c < args.distribution.end();
             ++c) {
            // Fill genotype tables
            for (i = 0; i < args.order - 2; ++i) {
                GenotypeTable<uint64_t>::combine(
                    i > 0 ? gts[i - 1] : args.dataset[c[0]],
                    args.dataset[c[i + 1]], gts[i]);
            }
            // Iterate over subsequent combinations
            for (i = c->back() + 1; i < (int)args.dataset.snps; ++i) {
                // If the block is full, compute all MI's
                if (j == B) {
                    for (k = 0; k < B; ++k) {
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
                    args.order > 2 ? gts.back() : args.dataset[c[0]],
                    args.dataset[i], cts[j++]);
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

    static void nonsegmented_search(Args &args)
    {
        int i, j, k;
        // Allocate genotype tables of size < target interaction order
        std::vector<GenotypeTable<uint64_t>> gts;
        gts.reserve(args.order - 2);
        for (auto o = 2; o < args.order; ++o) {
            gts.emplace_back(o, args.dataset[0].cases_words,
                             args.dataset[0].ctrls_words);
        }
        ContingencyTable<uint32_t> ct(args.order);
        Result<int, float> r;
        r.combination.resize(args.order);
        MutualInformation<float> mi(args.dataset.cases, args.dataset.ctrls);
        // For each combination assigned by the distribution
        for (auto c = args.distribution.begin(); c < args.distribution.end();
             ++c) {
            // Fill genotype tables
            for (i = 0; i < args.order - 2; ++i) {
                GenotypeTable<uint64_t>::combine(
                    i > 0 ? gts[i - 1] : args.dataset[c[0]],
                    args.dataset[c[i + 1]], gts[i]);
            }
#ifdef BENCHMARK
            args.combinations += args.dataset.snps - (c->back() + 1);
#endif
            // Iterate over subsequent combinations
            for (i = c->back() + 1; i < (int)args.dataset.snps; ++i) {
                memcpy(r.combination.data(), c->data(),
                       c->size() * sizeof(int));
                r.combination.back() = i;
                // Fill contingency table
                GenotypeTable<uint64_t>::combine_and_popcnt(
                    args.order > 2 ? gts.back() : args.dataset[c[0]],
                    args.dataset[i], ct);
                // Compute mutual information
                r.val = mi.compute(ct);
                args.maxarray.add(r);
            }
        }
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
