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
 * @file Arg_parser.cpp
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief Arg_parser class members implementation.
 */

#include <args.hxx>
#include <fiuncho/Definitions.h>
#include <fiuncho/utils/Arg_parser.h>
#include <fiuncho/utils/IOMpi.h>

std::istream &operator>>(std::istream &is,
                         std::pair<unsigned int, unsigned int> &ints) {
    is >> std::get<0>(ints);
    is.get();
    if (is.eof())
        throw args::ParseError("Pair missing second value");
    is >> std::get<1>(ints);
    return is;
}

Arg_parser::Arg_parser(int &argc, char **argv) : argc(argc), argv(argv) {}

Arg_parser::Arguments Arg_parser::get_arguments() {
    args::ArgumentParser parser(PROJECT_DESCRIPTION, PROJECT_HELP);
    args::Group required(parser, "Required arguments",
                         args::Group::Validators::All);
    args::Positional<std::string> r_tped(required, "tped-file",
                                         "path to TPED file");
    args::Positional<std::string> r_tfam(required, "tfam-file",
                                         "path to TFAM file");
    args::Positional<std::string> r_output(required, "output-file",
                                           "path to output file");
    args::Group cpu_opt(parser, "CPU runtime configuration",
                        args::Group::Validators::DontCare);
    args::ValueFlag<unsigned int> cpu_threads(
        cpu_opt, "thread-num", "number of threads to use per process",
        {'t', "threads"});
    args::Group verb(parser, "Verbosity level",
                     args::Group::Validators::AtMostOne);
    args::Flag verb_b(verb, "benchmarking", "print runtimes", {"benchmarking"});
    args::Flag verb_d(verb, "debug", "print debug information", {"debug"});
    args::ValueFlag<unsigned int> output_num(
        parser, "output-num", "number of triplets to print in the output file",
        {'n', "num-out"});
    args::ValueFlag<bool> use_mi(
        parser, "mutual-info",
        "use Mutual Information (the alternative is Information Gain,"
        " default = 1)",
        {"mi"}, true);
    // args::VersionFlag version(parser, "version", "output version information
    // and exit", {'V', "version"});
    args::HelpFlag help(parser, "help", "display this help and exit",
                        {'h', "help"});

    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        IOMpi::Instance().mprint<IOMpi::N>(parser.Help());
        throw Arg_parser::Finalize("hola mundo");
    } catch (const args::ParseError &e) {
        IOMpi::Instance().smprint<IOMpi::N>(std::cerr, e.what());
        IOMpi::Instance().smprint<IOMpi::N>(std::cerr, parser.Help());
        throw Arg_parser::Finalize("hola mundo");
    } catch (const args::ValidationError &e) {
        IOMpi::Instance().smprint<IOMpi::N>(std::cerr, e.what());
        IOMpi::Instance().smprint<IOMpi::N>(std::cerr, parser.Help());
        throw Arg_parser::Finalize("hola mundo");
    }

    Arguments arguments;
    arguments.tped = args::get(r_tped);
    arguments.tfam = args::get(r_tfam);
    arguments.output = args::get(r_output);
    arguments.benchmarking = verb_b;
    arguments.debug = verb_d;

    if (cpu_threads) {
        arguments.cpu_threads = args::get(cpu_threads);
    }
    if (output_num) {
        arguments.output_num = args::get(output_num);
    }
    if (use_mi) {
        arguments.use_mi = args::get(use_mi);
    }
    return arguments;
}