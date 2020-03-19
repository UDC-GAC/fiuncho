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
 * @file Engine.h
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief Abstract class Engine definition.
 */

#ifndef MPI3SNP_ENGINE_H
#define MPI3SNP_ENGINE_H

#include <vector>
#include "Statistics.h"
#include "Position.h"

class Engine {
public:
    class Error : public std::runtime_error {
    public:
        explicit Error(const std::string &message) : runtime_error(message) {};

        ~Error() override = default;
    };

    virtual void run(std::string tped, std::string tfam, std::vector<Position> &mutual_info, size_t num_outputs) = 0;
};

#endif //MPI3SNP_ENGINE_H
