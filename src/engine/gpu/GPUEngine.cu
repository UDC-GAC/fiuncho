/*
 * MPI3SNP ~ https://github.com/chponte/mpi3snp
 *
 * Copyright 2018 Christian Ponte
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
 * documentation files (the "Software"), to deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the
 * Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
 * WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * @file gpu/GPUEngine.cpp
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief GPUEngine class implementation.
 */

#include "GPUEngine.h"
#include "Distributor.h"
#include "MutualInformation.h"
#include <cstring>

GPUEngine::GPUEngine(unsigned int proc_num, unsigned int proc_id,
                     std::vector<std::pair<unsigned int, unsigned int>> gpu_map, bool use_mi, Statistics &statistics) :
        proc_num(proc_num),
        proc_id(proc_id),
        use_mi(use_mi),
        statistics(statistics) {
    statistics.Begin_timer("CUDA initialization time");
    cudaFree(nullptr);
    statistics.End_timer("CUDA initialization time");

    int avail_gpus = 0;
    if (cudaSuccess != cudaGetDeviceCount(&avail_gpus))
        throw CUDAError();
    if (avail_gpus == 0) {
        throw CUDAError("Could not find any CUDA-enabled GPU");
    }

    auto pos = std::find_if(gpu_map.begin(), gpu_map.end(),
                            [&proc_id](std::pair<unsigned int, unsigned int> item) { return item.first == proc_id; });
    gpu_id = pos == gpu_map.end() ? proc_id % avail_gpus : pos->second;

    cudaDeviceProp gpu_prop;
    if (cudaSuccess != cudaGetDeviceProperties(&gpu_prop, gpu_id))
        throw CUDAError();
    if (gpu_prop.major < 2 || !gpu_prop.canMapHostMemory) {
        throw CUDAError("GPU " + std::to_string(gpu_id) + " does not meet compute capabilities\n" +
                        "Name: " + gpu_prop.name + "\n" + "Compute capability: " +
                        std::to_string(gpu_prop.major) + "." + std::to_string(gpu_prop.minor));
    }
    if (cudaSuccess != cudaSetDevice(gpu_id))
        throw CUDAError();
}

void GPUEngine::run(std::string tped, std::string tfam, std::vector<Position> &output, size_t num_outputs) {
    statistics.Begin_timer("SNPs read time");
    Dataset *dataset;
    try {
        dataset = new Dataset(tped, tfam, Dataset::Transposed);
    } catch (const Dataset::Read_error &error) {
        throw Engine::Error(error.what());
    }
    statistics.End_timer("SNPs read time");

    statistics.Addi("SNP count", dataset->get_SNP_count());
    statistics.Addi("Number of cases", dataset->get_case_count());
    statistics.Addi("Number of controls", dataset->get_ctrl_count());

    Distributor<uint32_t, uint2> distributor(dataset->get_SNP_count(), proc_num);

    Algorithm<uint2> *search = new MutualInformation(use_mi, dataset->get_SNP_count(), dataset->get_case_count(), dataset->get_ctrl_count(),
                         dataset->get_cases(), dataset->get_ctrls());

    std::vector<uint2> pairs;
    distributor.get_pairs([](uint32_t x, uint32_t y) {
        uint2 p {x, y};
        return p;
    }, proc_id, pairs);

    std::string timer_label;
    timer_label += "GPU " + std::to_string(gpu_id) + " runtime";
    statistics.Begin_timer(timer_label);

    output.resize(num_outputs);
    long myTotalAnal = search->compute(pairs, num_outputs, &output.at(0));
    cudaDeviceSynchronize();
    statistics.Addl("GPU " + std::to_string(gpu_id) + " computations", myTotalAnal);

    delete search;
    delete dataset;

    statistics.End_timer(timer_label);
}