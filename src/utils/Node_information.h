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
 * @file Node_information.h
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief Abstract class definition. Node_information provides debug information about the nodes where MPI is running,
 * and the process allocation.
 */

#ifndef MPI3SNP_NODE_INFORMATION_H
#define MPI3SNP_NODE_INFORMATION_H

#include <vector>
#include <string>

class Node_information {
public:
    class Builder {
    public:
        static Node_information *get_information();
    };

    Node_information();

    virtual std::string hostname() const =0;

    virtual std::string cpu() const =0;

    virtual long memory() const =0;

    virtual std::string mpi_library_version() const =0;

    virtual std::vector<int> processes() const =0;

    virtual std::vector<std::string> gpus() const =0;

    virtual std::string to_string() const =0;

    static std::vector<Node_information *> gather(int process);

protected:
    virtual size_t to_byteblock(void **ptr) const =0;

    virtual void add_processes(std::vector<int> processes) =0;

    static Node_information *build_from_byteblock(const void *ptr);

    std::string hardware_id;

private:
    static constexpr const char *hardware_id_file = "/proc/net/netlink";
};


#endif //MPI3SNP_NODE_INFORMATION_H
