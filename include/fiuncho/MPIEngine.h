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
 * @file MPIEngine.h
 * @author Christian Ponte
 */

#ifndef FIUNCHO_MPIENGINE_H
#define FIUNCHO_MPIENGINE_H

#include <fiuncho/Search.h>
#include <fiuncho/utils/Result.h>
#include <limits>
#include <mpi.h>
#include <sstream>
#include <string>
#include <vector>

#ifdef BENCHMARK
#include <iostream>
#endif

/**
 * Epistasis search engine, implementing a distributed algorithm using MPI. Each
 * node available in the MPI context is used to explore SNP combinations in
 * parallel.
 */

class MPIEngine
{

    const int mpi_size;
    const int mpi_rank;

    int get_mpi_size()
    {
        int size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        return size;
    }

    int get_mpi_rank()
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        return rank;
    }

    static std::string
    serialize_results(const std::vector<Result<int, float>> &results)
    {
        std::stringstream oss;
        for (auto r = results.begin(); r != results.end(); r++) {
            Result<int, float>::serialize(oss, *r);
        }
        std::string s = oss.str();
        return s;
    }

    static void deserialize_results(const std::string &s,
                                    std::vector<Result<int, float>> &v)
    {
        std::stringstream iss(s);
        Result<int, float> buffer;
        while ((size_t)iss.tellg() < s.size()) {
            Result<int, float>::deserialize(iss, buffer);
            v.push_back(buffer);
        }
    }

  public:
    /**
     * @name Constructors
     */
    //@{

    /**
     * Create an MPIEngine object. The constructor calls MPI routines, and thus
     * it is mandatory to call the constructor after the MPI environment has
     * been initialized with the `MPI_Init` function.
     */

    MPIEngine() : mpi_size(get_mpi_size()), mpi_rank(get_mpi_rank()) {}

    //@}

    /**
     * @name Methods
     */
    //@{

    /**
     * Run the epistasis search on the different MPI processes. Each process
     * will, in turn, call Search::run to exploit the resources available to
     * that process. The returned vector will only be available to process 0.
     *
     * @return Vector of Result's sorted in descending order by their
     * MutualInformation value
     * @param tped Path to the tped data file
     * @param tfam Path to the tfam data file
     * @param order Order of the epistatic interactions to locate
     * @param outputs Number of results to include in the output vector
     * @param args Arguments to the Search class
     * @tparam T Search class to use in the epistasis search
     * @tparam Args Argument types of the Search class constructor. This
     * template parameter should be automatically deduced by the compiler and
     * its explicit use is discouraged
     */

    template <typename T, typename... Args>
    std::vector<Result<int, float>>
    run(const std::string &tped, const std::string &tfam,
        const unsigned int order, const unsigned int outputs, Args &&...args)
    {
        std::vector<Result<int, float>> local_results, global_results;
#ifdef BENCHMARK
        double function_time, dataset_time;
        function_time = MPI_Wtime();
        dataset_time = MPI_Wtime();
#endif
#ifdef ALIGN
        const auto dataset = Dataset<uint64_t>::read<ALIGN>(tped, tfam);
#else
        const auto dataset = Dataset<uint64_t>::read(tped, tfam);
#endif
        // Check Dataset size to avoid int overflow
        if (dataset.snps > (size_t)std::numeric_limits<int>::max()) {
            throw std::runtime_error(
                "Input data limit exceeded: Dataset contains more than " +
                std::to_string(std::numeric_limits<int>::max()) + " SNPs");
        }
#ifdef BENCHMARK
        dataset_time = MPI_Wtime() - dataset_time;
        std::cout << "Read " << dataset.snps << " SNPs from "
                  << dataset.cases + dataset.ctrls << " individuals ("
                  << dataset.cases << " cases, " << dataset.ctrls
                  << " controls) in " << dataset_time << " seconds\n";
#endif
        const Distribution<int> distribution(dataset.snps, order - 1, mpi_size,
                                             mpi_rank);
        Search *search = new T(std::forward<Args>(args)...);
        local_results = search->run(dataset, order, distribution, outputs);
        delete search;
        // Serialize the results
        const std::string serialized_results = serialize_results(local_results);
        local_results.clear();
        const int nbytes = serialized_results.size();
        // If process is rank 0
        if (mpi_rank == 0) {
            // Gather the number of results per process
            std::vector<int> recv_counts(mpi_size);
            MPI_Gather(&nbytes, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, 0,
                       MPI_COMM_WORLD);
            std::vector<int> recv_displs(mpi_size);
            recv_displs[0] = 0;
            for (auto i = 1; i < mpi_size; i++) {
                recv_displs[i] = recv_displs[i - 1] + recv_counts[i - 1];
            }
            // Gather the actual results
            std::string buffer(
                recv_displs[mpi_size - 1] + recv_counts[mpi_size - 1], '\0');
            MPI_Gatherv(serialized_results.data(), serialized_results.size(),
                        MPI_BYTE, (void *)buffer.data(), recv_counts.data(),
                        recv_displs.data(), MPI_BYTE, 0, MPI_COMM_WORLD);
            // Deserialize the results
            deserialize_results(buffer, global_results);
            // Sort the result
            std::sort(global_results.rbegin(), global_results.rend());
            if (global_results.size() > outputs) {
                global_results.resize(outputs);
            }
        } else {
            // Send the number of results to process with rank 0
            MPI_Gather(&nbytes, 1, MPI_INT, nullptr, 1, MPI_INT, 0,
                       MPI_COMM_WORLD);
            // Send the actual serialized results to process with rank 0
            MPI_Gatherv(serialized_results.data(), serialized_results.size(),
                        MPI_BYTE, nullptr, nullptr, nullptr, MPI_BYTE, 0,
                        MPI_COMM_WORLD);
        }

#ifdef BENCHMARK
        function_time = MPI_Wtime() - function_time;
        std::cout << "Total elapsed time: " << function_time << '\n';
#endif

        return global_results;
    }

    //@}
};

#endif
