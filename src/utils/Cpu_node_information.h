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
 * @file Cpu_node_information.h
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief Declaration of the abstract class Node_information, specific for CPU architectures. Provides the host name,
 * cpu name, memory size, specific mpi library implementation and the list of processes hosted in that node.
 */

#ifndef MPI3SNP_CPU_NODE_INFORMATION_H
#define MPI3SNP_CPU_NODE_INFORMATION_H

#include "Node_information.h"

class Cpu_node_information : public Node_information {
public:
    Cpu_node_information();

    explicit Cpu_node_information(const void *ptr);

    std::string hostname() const override;

    std::string cpu() const override;

    long memory() const override;

    std::string mpi_library_version() const override;

    std::vector<int> processes() const override;

    std::vector<std::string> gpus() const override;

    std::string to_string() const override;

protected:
    size_t to_byteblock(void **ptr) const override;

    void add_processes(std::vector<int> processes) override;

private:
    static constexpr const char *cpu_file = "/proc/cpuinfo";

    static constexpr const char *memory_file = "/proc/meminfo";

    static std::string get_cpu_info();

    static long get_physical_memory();

    std::string hostname_str;
    std::string cpu_str;
    long memory_size;
    std::string mpi_library;
    std::vector<int> process_list;
};


#endif //MPI3SNP_CPU_NODE_INFORMATION_H
