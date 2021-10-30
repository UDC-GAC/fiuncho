/**
 * @file main.cpp
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief Main program. Initializes the MPI environment, parses the command-line
 * arguments, prints debug information if necessary, instantiates the Search and
 * terminates the execution.
 */

#include <filesystem>
#include <fiuncho/MPIEngine.h>
#include <fiuncho/ThreadedSearch.h>
#include <fstream>
#include <iostream>
#include <tclap/CmdLine.h>
#include <unistd.h>

typedef struct {
    std::vector<std::string> inputs;
    std::string output;
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
    TCLAP::UnlabeledMultiArg<std::string> files(
        "files", "input and output files (output goes last)", true,
        "list of string");
    cmd.add(files);
    // Read
    Arguments args;
    cmd.parse(argc, argv);
    args.inputs = files.getValue();
    args.output = args.inputs.back();
    args.inputs.pop_back();
    // Check that all input files are readable
    std::for_each(args.inputs.begin(), args.inputs.end(),
                  [](const std::string &path) {
                      auto absolute = std::filesystem::absolute(path);
                      if (access(absolute.string().c_str(), R_OK) != 0)
                          throw std::runtime_error("Cannot read file \"" +
                                                   absolute.string() + "\"");
                  });
    // Check that the output file is writeable
    [](const std::string &path) {
        auto absolute = std::filesystem::absolute(path);
        auto status = std::filesystem::status(absolute);
        // If the file exists, check if we have write perm
        if (status.type() == std::filesystem::file_type::regular) {
            if (access(absolute.string().c_str(), W_OK) != 0)
                throw std::runtime_error("File \"" + absolute.string() +
                                         "\" is not writeable");
        } else { // File does not exist
            auto parent = std::filesystem::path(absolute).parent_path();
            if (access(parent.string().c_str(), W_OK) != 0)
                throw std::runtime_error("Cannot write into directory \"" +
                                         parent.string() + "\"");
        }
    }(args.output);
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
        auto results = engine.run<ThreadedSearch>(args.inputs, args.order,
                                                  args.noutputs, args.threads);
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
