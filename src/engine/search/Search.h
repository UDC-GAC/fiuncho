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
 * @file Search.h
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief Search class definition.
 */

#ifndef MPI3SNP_SEARCH_H
#define MPI3SNP_SEARCH_H

#include <vector>
#include "Arg_parser.h"
#include "Engine.h"

class Search {
public:
    class Builder {
    public:
        static Search *build_from_args(Arg_parser::Arguments arguments, Statistics &statistics);

        Builder() = delete;
    };

    // Execute the epistasis search
    void execute();

private:
    Search() = default;

    Engine *engine;
    int proc_id, num_proc;
    std::string tped_file, tfam_file, out_file;
    unsigned int num_outputs;
};

#endif //MPI3SNP_SEARCH_H
