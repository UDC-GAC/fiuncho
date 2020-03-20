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
 * @file Arg_parser.h
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief Argument parser class. Takes the command line arguments as the input,
 * and returns a structure containing the configuration to be used.
 */

#ifndef FIUNCHO_ARG_PARSER_H
#define FIUNCHO_ARG_PARSER_H

#include <map>
#include <thread>
#include <vector>

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
        explicit Finalize(const std::string &message)
            : runtime_error(message){};

        ~Finalize() override = default;
    };

    Arg_parser(int &argc, char **argv);

    Arguments get_arguments();

  private:
    const int argc;
    char **argv;
};

#endif
