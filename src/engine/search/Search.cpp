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
 * @file Search.cpp
 * @author Christian Ponte
 * @date 1 October 2018
 *
 * @brief Search class members implementation.
 */

#include "Search.h"
#include "IOMpi.h"
#include "Dataset.h"
#include "Definitions.h"
#include <fstream>
#include <mpi.h>

#if TARGET_ARCH == CPU

#include "CPUEngine.h"

#else

#include "GPUEngine.h"
#include "CUDAError.h"

#endif

Search *Search::Builder::build_from_args(Arg_parser::Arguments arguments, Statistics &statistics) {
    auto *search = new Search();
    MPI_Comm_rank(MPI_COMM_WORLD, &search->proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &search->num_proc);
    search->tped_file = arguments.tped;
    search->tfam_file = arguments.tfam;
    search->out_file = arguments.output;
    search->num_outputs = arguments.output_num;

#if TARGET_ARCH == CPU
    search->engine = new CPUEngine(search->num_proc, search->proc_id, arguments.cpu_threads, arguments.use_mi,
                                   statistics);
#else
    try {
        search->engine = new GPUEngine(search->num_proc, search->proc_id, arguments.gpu_map, arguments.use_mi, statistics);
    } catch (const Engine::Error &e) {
        IOMpi::Instance().smprint<IOMpi::E>(std::cerr, std::string(e.what()) + "\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
        return nullptr;
    }
#endif
    return search;
}

void Search::execute() {
    std::vector<Position> mutual_info, result;

    try {
        engine->run(tped_file, tfam_file, mutual_info, num_outputs);
    } catch (const Engine::Error &e) {
        IOMpi::Instance().smprint<IOMpi::E>(std::cerr, std::string(e.what()) + "\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
        return;
    }

    // Gather the results on the master process
    if (proc_id == 0) {
        result.resize(num_outputs * num_proc);
    }

    MPI_Gather(&mutual_info.front(), num_outputs * sizeof(Position), MPI_BYTE,
               &result.front(), num_outputs * sizeof(Position), MPI_BYTE, 0, MPI_COMM_WORLD);

    if (proc_id == 0) {
        // Sort the result
        std::sort(result.begin(), result.end(), [](Position a, Position b) { return b < a; });

        // Write results to the output file
        std::ofstream of(out_file.c_str(), std::ios::out);
        for (unsigned int i = 0; i < num_outputs; i++) {
            of << result[i].to_string() << '\n';
        }
        of.close();
    }

    IOMpi::Instance().cprint<IOMpi::D>("3-SNP analysis finalized\n");
}