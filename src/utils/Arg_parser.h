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
 * @file Arg_parser.h
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief Argument parser class. Takes the command line arguments as the input, and returns a structure containing the
 * configuration to be used.
 */

#ifndef MPI3SNP_ARG_PARSER_H
#define MPI3SNP_ARG_PARSER_H

#include <map>
#include <vector>
#include <thread>

class Arg_parser {
public:
    struct Arguments {
        std::string tped;
        std::string tfam;
        std::string output;
        std::vector<std::pair<unsigned int, unsigned int>> gpu_map;
        unsigned int cpu_threads = std::thread::hardware_concurrency();
        bool use_mi = true;
        unsigned int output_num = 10;
        bool debug;
        bool benchmarking;
    };

    class Finalize : public std::runtime_error {
    public:
        explicit Finalize(const std::string &message) : runtime_error(message) {};

        ~Finalize() override = default;
    };

    Arg_parser(int &argc, char **argv);

    Arguments get_arguments();

private:
    const int argc;
    char **argv;
};

#endif //MPI3SNP_ARG_PARSER_H
