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
 * @file gpu/Gpu_node_information.h
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief GPU implementation of the abstract class Node_information, extending the CPU implementation provided in
 * Cpu_node_information. This class adds information about the GPUs present in the node.
 */

#ifndef MPI3SNP_GPU_NODE_INFORMATION_H
#define MPI3SNP_GPU_NODE_INFORMATION_H

#include "../Cpu_node_information.h"

class Gpu_node_information : public Cpu_node_information {
public:
    Gpu_node_information();

    explicit Gpu_node_information(const void *ptr);

    std::vector<std::string> gpus() const override;

    std::string to_string() const override;

protected:
    size_t to_byteblock(void **ptr) const override;

private:
    std::vector<std::string> gpu_list;
};

#endif //MPI3SNP_GPU_NODE_INFORMATION_H
