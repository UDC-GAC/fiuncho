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
#include <cstdint>
#include <fiuncho/GenotypeTable.h>
#include <fiuncho/dataset/Dataset.h>
#include <iostream>
#include <pthread.h>
#include <thread>
#include <time.h>

/*
 *  This benchmark programs calculates bit tables of a specific order following
 *  this process:
 *      1. Spawn as many threads as indicated
 *      2. Set processor affinity as indicated
 *      3. Initialize previous tables and such, and synchronize threads at the
 *         end of the initialization
 *      4. Warm up the CPU core by running an additional 10% of the total
 *         iterations
 *      5. Measure bit table computation time
 *      6. Print elapsed time on each thread
 *
 *  Program arguments:
 *      1: Comma-separated list of hardware threads to run on.
 *      2: Order of the tables to benchmark
 *      3: How many times the main loop is repeated
 *      4: Path to the TPED input file
 *      5: Path to the TFAM input file
 */

int repetitions;

unsigned short thread_count;
pthread_barrier_t barrier;

void bench(const std::string tped, const std::string tfam,
           const unsigned short order, double &elapsed_time, const int affinity)
{
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(affinity, &cpuset);
    int rc = pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
    if (rc != 0) {
        std::cerr << "Error calling pthread_setaffinity_np: " << rc << "\n";
    }

#ifdef ALIGN
    const auto dataset = Dataset<uint64_t>::read<ALIGN>(tped, tfam);
#else
    const auto dataset = Dataset<uint64_t>::read(tped, tfam);
#endif

    const size_t snp_count = dataset.snps;
    GenotypeTable<uint64_t> table(order, dataset[0].cases_words,
                                  dataset[0].ctrls_words);
    struct timespec start, end;

    if (order == 2) {
        // Main compute loop for order == 2
        pthread_barrier_wait(&barrier);
        // Warmup CPU adding an extra 10% of iterations before measuring time
        for (auto reps = 0; reps < repetitions / 10; reps++) {
            for (size_t snp = 1; snp < snp_count; snp++) {
                GenotypeTable<uint64_t>::combine(dataset[0], dataset[snp],
                                                 table);
            }
        }
        // Measure time
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start);
        for (auto reps = 0; reps < repetitions; reps++) {
            for (size_t snp = 1; snp < snp_count; snp++) {
                GenotypeTable<uint64_t>::combine(dataset[0], dataset[snp],
                                                 table);
            }
        }
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end);
        elapsed_time = end.tv_sec + end.tv_nsec * 1E-9 - start.tv_sec -
                       start.tv_nsec * 1E-9;
    } else {
        // Initialize previous bit tables outside of the main loop
        std::vector<GenotypeTable<uint64_t>> prev_tables;
        prev_tables.emplace_back(2, dataset[0].cases_words,
                                 dataset[0].ctrls_words);
        GenotypeTable<uint64_t>::combine(dataset[0], dataset[1],
                                         prev_tables[0]);
        for (auto o = 3; o < order; o++) {
            prev_tables.emplace_back(o, dataset[0].cases_words,
                                     dataset[0].ctrls_words);
            GenotypeTable<uint64_t>::combine(prev_tables[o - 3], dataset[o - 1],
                                             prev_tables[o - 2]);
        }
        const GenotypeTable<uint64_t> &last = prev_tables[order - 3];
        // Main compute loop for order > 2
        pthread_barrier_wait(&barrier);
        // Warmup CPU adding an extra 10% of iterations before measuring time
        for (auto reps = 0; reps < repetitions / 10; reps++) {
            for (size_t snp = order - 1; snp < snp_count; snp++) {
                GenotypeTable<uint64_t>::combine(last, dataset[snp], table);
            }
        }
        // Measure time
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start);
        for (auto reps = 0; reps < repetitions; reps++) {
            for (size_t snp = order - 1; snp < snp_count; snp++) {
                GenotypeTable<uint64_t>::combine(last, dataset[snp], table);
            }
        }
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end);
        elapsed_time = end.tv_sec + end.tv_nsec * 1E-9 - start.tv_sec -
                       start.tv_nsec * 1E-9;
    }
}

int main(int argc, char *argv[])
{
    if (argc != 6) {
        std::cout << argv[0] << " <THREADS> <ORDER> <REPETITIONS> <TPED> <TFAM>"
                  << std::endl;
        return 0;
    }

    // Initialization
    // Arguments
    std::vector<int> affinity = split_into_ints(std::string(argv[1]), ',');
    thread_count = affinity.size();
    const unsigned short order = atoi(argv[2]);
    repetitions = atoi(argv[3]);
    const std::string tped = argv[4], tfam = argv[5];
    // Variables
    pthread_barrier_init(&barrier, NULL, thread_count);
    std::vector<std::thread> threads;
    std::vector<double> times;
    times.resize(thread_count);

    // Spawn thread_count - 1 threads
    for (size_t i = 1; i < affinity.size(); i++) {
        threads.emplace_back(bench, tped, tfam, order, std::ref(times[i]),
                             affinity[i]);
    }
    // Also use current thread
    bench(tped, tfam, order, times[0], affinity[0]);

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
