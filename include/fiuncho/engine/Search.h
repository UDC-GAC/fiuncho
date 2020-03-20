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
 * @file Search.h
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief Search class definition.
 */

#include <fiuncho/engine/Engine.h>
#include <fiuncho/utils/Arg_parser.h>
#include <vector>

#ifndef FIUNCHO_SEARCH_H
#define FIUNCHO_SEARCH_H

class Search {
  public:
    class Builder {
      public:
        static Search *build_from_args(Arg_parser::Arguments arguments,
                                       Statistics &statistics);

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

#endif
