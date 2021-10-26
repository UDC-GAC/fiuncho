/**
 * @file gtable.cpp
 * @author Christian (christian.ponte@udc.es)
 * @brief Benchmarking program for the genotype table computation function
 * @date 20/09/2021
 */

#include "utils.h"
#include <fiuncho/GenotypeTable.h>
#include <fiuncho/dataset/Dataset.h>
#include <iostream>

/**
 * @brief Genotype table calculation loop. The function repeats, for a given
 * number of iterations, the computation of the genotype tables for {SNP_0, ...,
 * SNP_o-2, SNP_i}, with o being the order of the genotype tables and i in [o-1,
 * dataset.snps].
 *
 * @param repetitions Number of iterations to loop
 * @param dataset Input data
 * @param order Order of the genotype tables computed
 * @param prev Genotype table of combination {SNP_0, ..., SNP_o-2}
 */

void loop_gtable(const int repetitions, const Dataset<uint64_t> &dataset,
                 const unsigned short order,
                 const GenotypeTable<uint64_t> &prev)
{
    const size_t snp_count = dataset.snps;
    GenotypeTable<uint64_t> gtable(order, dataset[0].cases_words,
                                   dataset[0].ctrls_words);
    for (auto reps = 0; reps < repetitions; reps++) {
        for (size_t snp = order - 1; snp < snp_count; snp++) {
            GenotypeTable<uint64_t>::combine(prev, dataset[snp], gtable);
        }
    }
}

/**
 * @brief Main benchmarking function. It precomputes the last genotype table
 * corresponding to the combination {SNP_0, ..., SNP_o-2} and measures the
 * elapsed time during the calculation of the genotype table of the combination
 * {SNP_0, ..., SNP_o-2, SNP_i}, with o being the order of the genotype tables
 * and i in [o-1, dataset.snps]. The function creates as many threads as CPU
 * cores indicated in the affinity, and pins them to the specified allocation.
 *
 * @param tped TPED file path
 * @param tfam TFAM file path
 * @param order Order of the genotype tables computed by the benchmark
 * @param affinity List of CPU cores to be used
 * @param repetitions Number of repetitions to loop during the measurement
 * @return A vector containing the elapsed times for each thread
 */

std::vector<double> bench_gtable(const std::string tped, const std::string tfam,
                                 const unsigned short order,
                                 const std::vector<int> &affinity,
                                 const int repetitions)
{
    const auto dataset = Dataset<uint64_t>::read(tped, tfam);
    const size_t cases_words = dataset[0].cases_words,
                 ctrls_words = dataset[1].ctrls_words;
    if (order == 2) {
        return multithread_pinned_time(affinity, loop_gtable, repetitions,
                                       dataset, order, dataset[0]);
    } else {
        // Initialize previous bit tables outside of the main loop
        std::vector<GenotypeTable<uint64_t>> gtables;
        gtables.emplace_back(2, cases_words, ctrls_words);
        GenotypeTable<uint64_t>::combine(dataset[0], dataset[1], gtables[0]);
        for (auto o = 3; o < order; o++) {
            gtables.emplace_back(o, cases_words, ctrls_words);
            GenotypeTable<uint64_t>::combine(gtables[o - 3], dataset[o - 1],
                                             gtables[o - 2]);
        }
        return multithread_pinned_time(affinity, loop_gtable, repetitions,
                                       dataset, order, gtables.back());
    }
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
    auto times = bench_gtable(tped, tfam, order, affinity, repetitions);
    // Print times
    for (size_t i = 0; i < times.size() - 1; i++) {
        std::cout << times[i] << ',';
    }
    std::cout << times.back() << '\n';

    return 0;
}
