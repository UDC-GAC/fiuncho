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
 * @file gpu/Gpu_node_information.cpp
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief Gpu_node_information class members implementation.
 */

#include "Gpu_node_information.h"
#include <algorithm>
#include <cstring>
#include <cuda_runtime.h>
#include <iostream>

Gpu_node_information::Gpu_node_information() : Cpu_node_information() {
    int avail_gpus = 0;
    if (cudaSuccess != cudaGetDeviceCount(&avail_gpus)) {
        return;
    }

    cudaDeviceProp gpu_prop;
    for (int gid = 0; gid < avail_gpus; gid++) {
        if (cudaSuccess != cudaGetDeviceProperties(&gpu_prop, gid)) {
            continue;
        }
        gpu_list.push_back(std::to_string(gid) + ": " + gpu_prop.name);
    }
}

Gpu_node_information::Gpu_node_information(const void *ptr) : Cpu_node_information(((char *) ptr) + sizeof(size_t)) {
    auto *block = (char *) ptr;
    size_t size, offset = 0;
    memcpy(&offset, block, sizeof(size_t));

    // GPU information vector
    memcpy(&size, block + offset, sizeof(size_t));
    offset += sizeof(size_t);
    gpu_list.resize(size);
    for (auto &g : gpu_list) {
        memcpy(&size, block + offset, sizeof(size_t));
        offset += sizeof(size_t);
        char gpu_info[size];
        memcpy(gpu_info, block + offset, size);
        offset += size;
        g = std::string(gpu_info, size);
    }
}

std::vector<std::string> Gpu_node_information::gpus() const {
    return gpu_list;
}

std::string Gpu_node_information::to_string() const {
    std::string output = Cpu_node_information::to_string();
    output += "GPUs: ";
    for (const auto &g : gpu_list){
        output += "[" + g + "] ";
    }
    output += "\n";
    return output;
}

size_t Gpu_node_information::to_byteblock(void **ptr) const {
    void *cpu_info;
    size_t cpu_block_size = Cpu_node_information::to_byteblock(&cpu_info);
    // Size = number of strings + (length of string + string) for each string in the vector
    size_t gpu_vector_size = sizeof(size_t);
    std::for_each(gpu_list.begin(), gpu_list.end(),
                  [&gpu_vector_size](auto str) { gpu_vector_size += sizeof(size_t) + str.length(); });
    // Memory allocation
    *ptr = new char[sizeof(size_t) + cpu_block_size + gpu_vector_size];
    auto *buffer = (char *) *ptr;
    size_t size, offset = 0;

    // Length of CPU information byte block
    size = sizeof(size_t) + cpu_block_size;
    memcpy(buffer + offset, &size, sizeof(size_t));
    offset += sizeof(size_t);

    // CPU information byte block
    memcpy(buffer + offset, cpu_info, cpu_block_size);
    offset += cpu_block_size;

    // GPU information vector
    size = gpu_list.size();
    memcpy(buffer + offset, &size, sizeof(size_t));
    offset += sizeof(size_t);
    for (const auto &g : gpu_list) {
        size = g.length();
        memcpy(buffer + offset, &size, sizeof(size_t));
        offset += sizeof(size_t);
        memcpy(buffer + offset, &g[0], size);
        offset += size;
    }
    return offset;
}
