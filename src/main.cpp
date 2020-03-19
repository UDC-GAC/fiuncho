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
 * @file main.cpp
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief Main program. Initializes the MPI environment, parses the command-line arguments, prints debug information if
 * necessary, instantiates the Search and terminates the execution.
 */

#include <mpi.h>
#include "IOMpi.h"
#include "Search.h"
#include "Node_information.h"

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    double stime, etime;

    Arg_parser::Arguments arguments;
    try {
        arguments = Arg_parser(argc, argv).get_arguments();
    } catch (const Arg_parser::Finalize &finalize) {
        MPI_Finalize();
        return 0;
    }

    // Set adequate information printing level
    if (arguments.benchmarking) {
        IOMpi::Instance().Set_print_level(IOMpi::B);
    } else if (arguments.debug) {
        IOMpi::Instance().Set_print_level(IOMpi::D);

        // Print node information
        auto info_list = Node_information::gather(0);
        IOMpi::Instance().mprint<IOMpi::D>("*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*\n");
        for (int i = 0; i < info_list.size(); i++) {
            IOMpi::Instance().mprint<IOMpi::D>(info_list[i]->to_string());
            if (i == info_list.size() - 1) {
                IOMpi::Instance().mprint<IOMpi::D>("*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*\n\n");
            } else {
                IOMpi::Instance().mprint<IOMpi::D>("---------------------------------------------\n");
            }
            delete info_list[i];
        }
        info_list.clear();
    }

    /*get the startup time*/
    stime = MPI_Wtime();

    // Execute search
    Statistics statistics;
    Search *search = Search::Builder::build_from_args(arguments, statistics);
    if (search == nullptr) {
        MPI_Finalize();
        return 0;
    }
    search->execute();

    /*get the ending time*/
    etime = MPI_Wtime();

    delete search;
    // Print runtime statistics to stdout
    IOMpi::Instance().cprint<IOMpi::B>(statistics.To_string());
    IOMpi::Instance().mprint<IOMpi::B>("Overall time: " + std::to_string(etime - stime) + " seconds\n");

    MPI_Finalize();
    return 0;
}
