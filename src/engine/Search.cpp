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
 * @file Search.cpp
 * @author Christian Ponte
 * @date 1 October 2018
 *
 * @brief Search class members implementation.
 */

#include <fiuncho/Definitions.h>
#include <fiuncho/dataset/Dataset.h>
#include <fiuncho/engine/Search.h>
#include <fiuncho/engine/cpu/CPUEngine.h>
#include <fiuncho/utils/IOMpi.h>
#include <fstream>
#include <mpi.h>

Search *Search::Builder::build_from_args(Arg_parser::Arguments arguments,
                                         Statistics &statistics) {
    auto *search = new Search();
    MPI_Comm_rank(MPI_COMM_WORLD, &search->proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &search->num_proc);
    search->tped_file = arguments.tped;
    search->tfam_file = arguments.tfam;
    search->out_file = arguments.output;
    search->num_outputs = arguments.output_num;
    search->engine =
        new CPUEngine(search->num_proc, search->proc_id, arguments.cpu_threads,
                      arguments.use_mi, statistics);
    return search;
}

void Search::execute() {
    std::vector<Position> mutual_info, result;

    try {
        engine->run(tped_file, tfam_file, mutual_info, num_outputs);
    } catch (const Engine::Error &e) {
        IOMpi::Instance().smprint<IOMpi::E>(std::cerr,
                                            std::string(e.what()) + "\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
        return;
    }

    // Gather the results on the master process
    if (proc_id == 0) {
        result.resize(num_outputs * num_proc);
    }

    MPI_Gather(&mutual_info.front(), num_outputs * sizeof(Position), MPI_BYTE,
               &result.front(), num_outputs * sizeof(Position), MPI_BYTE, 0,
               MPI_COMM_WORLD);

    if (proc_id == 0) {
        // Sort the result
        std::sort(result.begin(), result.end(),
                  [](Position a, Position b) { return b < a; });

        // Write results to the output file
        std::ofstream of(out_file.c_str(), std::ios::out);
        for (unsigned int i = 0; i < num_outputs; i++) {
            of << result[i].to_string() << '\n';
        }
        of.close();
    }

    IOMpi::Instance().cprint<IOMpi::D>("3-SNP analysis finalized\n");
}