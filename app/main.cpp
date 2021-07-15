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
 * @file main.cpp
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief Main program. Initializes the MPI environment, parses the command-line
 * arguments, prints debug information if necessary, instantiates the Search and
 * terminates the execution.
 */

#include <fiuncho/MPIEngine.h>
#include <fiuncho/ThreadedSearch.h>
#include <fstream>
#include <iostream>
#include <linux/limits.h>
#include <tclap/CmdLine.h>
#include <unistd.h>

typedef struct {
    std::string tped, tfam, output;
    short order, threads;
    unsigned int noutputs;
} Arguments;

Arguments read_arguments(int argc, char **argv)
{
    // Create TCLAP CmdLine
    TCLAP::CmdLine cmd(
        "Full documentation available at https://fiuncho.readthedocs.io/", ' ',
        FIUNCHO_VERSION);
    class : public TCLAP::StdOutput
    {
      public:
        virtual void version(TCLAP::CmdLineInterface &c)
        {
            std::cout << "Fiuncho"
                      << " " << c.getVersion() << " (" << FIUNCHO_COMMIT_HASH
                      << ")\nBuilt with " << COMPILER_NAME << " ("
                      << COMPILER_VERSION << ") using flags: " << COMPILER_FLAGS
                      << "\n\n";
            std::cout
                << "Copyright (C) 2021 Christian Ponte, Universidade da "
                   "Coruña\n"
                   "This program comes with ABSOLUTELY NO WARRANTY. This is "
                   "free\n"
                   "software, and you are welcome to redistribute it under\n"
                   "certain conditions. License available at:\n"
                   "https://fiuncho.readthedocs.io/en/latest/license.html"
                << std::endl;
        }

        virtual void failure(TCLAP::CmdLineInterface &, TCLAP::ArgException &e)
        {
            throw e;
        }
    } cmd_output;
    cmd.setOutput(&cmd_output);
    class : public TCLAP::Constraint<short>
    {
        bool check(const short &order) const { return order > 1; }

        std::string shortID() const { return "integer"; }

        std::string description() const
        {
            return "order is equal or greater than 2";
        }
    } order_constraint;
    TCLAP::ValueArg<short> order("o", "order",
                                 "Order of the interactions to explore.", true,
                                 0, &order_constraint);
    cmd.add(order);
    class : public TCLAP::Constraint<short>
    {
        bool check(const short &threads) const { return threads > 0; }

        std::string shortID() const { return "integer"; }

        std::string description() const { return "threads is greater than 0"; }
    } threads_constraint;
    TCLAP::ValueArg<short> threads("t", "threads",
                                   "Number of threads to use. By default it "
                                   "uses all cores available to the process.",
                                   false, 0, &threads_constraint);
    cmd.add(threads);
    class : public TCLAP::Constraint<int>
    {
        bool check(const int &noutputs) const { return noutputs > 0; }

        std::string shortID() const { return "integer"; }

        std::string description() const { return "noutputs is greater than 0"; }
    } noutputs_constraint;
    TCLAP::ValueArg<int> noutputs("n", "noutputs",
                                  "Number of combinations to output. By "
                                  "default, it outputs 10 combinations.",
                                  false, 10, &noutputs_constraint);
    cmd.add(noutputs);
    class : public TCLAP::Constraint<std::string>
    {
        bool check(const std::string &path) const
        {
            return access(path.c_str(), R_OK) == 0;
        }

        std::string shortID() const { return "path"; }

        std::string description() const
        {
            return "path points to a readable file";
        }
    } infile_constraint;
    TCLAP::UnlabeledValueArg<std::string> tped(
        "tped", "Path to the input tped data file.", true, "",
        &infile_constraint);
    cmd.add(tped);
    TCLAP::UnlabeledValueArg<std::string> tfam(
        "tfam", "Tath to the input tfam data file.", true, "",
        &infile_constraint);
    cmd.add(tfam);
    class : public TCLAP::Constraint<std::string>
    {
        bool check(const std::string &path) const
        {
            if (access(path.c_str(), F_OK) == 0) { // File exists
                return access(path.c_str(), W_OK) == 0;
            } else { // File does not exist
                auto pos = path.find_last_of("/");
                // If there are no forward slashes in the path, the path points
                // to a file in the current directory
                if (pos == std::string::npos) {
                    char current_path[PATH_MAX + 1];
                    if (getcwd(current_path, PATH_MAX + 1) == nullptr) {
                        return false;
                    }
                    return access(current_path, W_OK) == 0;
                } else {
                    std::string parent_dir = path.substr(0, pos);
                    return access(parent_dir.c_str(), W_OK) == 0;
                }
            }
        }

        std::string shortID() const { return "path"; }

        std::string description() const
        {
            return "path points to a writeable location";
        }
    } outfile_constraint;
    TCLAP::UnlabeledValueArg<std::string> output(
        "output", "Path to the output file.", true, "", &outfile_constraint);
    cmd.add(output);
    // Read
    Arguments args;
    cmd.parse(argc, argv);
    args.tped = tped.getValue();
    args.tfam = tfam.getValue();
    args.output = output.getValue();
    args.order = order.getValue();
    args.threads = threads.getValue();
    args.noutputs = noutputs.getValue();
    return args;
}

int main(int argc, char **argv)
{
    int rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    try {
        // Read arguments
        auto args = read_arguments(argc, argv);
        // Execute search
        MPIEngine engine;
        auto results = engine.run<ThreadedSearch>(
            args.tped, args.tfam, args.order, args.noutputs, args.threads);
        if (rank == 0) {
            // Write results to the output file
            std::ofstream of(args.output, std::ios::out);
            for (auto r : results) {
                of << r.str() << '\n';
            }
            of.close();
        }
    } catch (const TCLAP::ArgException &e) {
        std::cerr << e.error() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    } catch (const std::runtime_error &e) {
        std::cerr << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Finalize();
    return 0;
}
