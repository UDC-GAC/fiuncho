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
 * @file Engine.h
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief Abstract class Engine definition.
 */

#ifndef FIUNCHO_ENGINE_H
#define FIUNCHO_ENGINE_H

#include <fiuncho/engine/Position.h>
#include <fiuncho/utils/Statistics.h>
#include <vector>

class Engine {
  public:
    class Error : public std::runtime_error {
      public:
        explicit Error(const std::string &message) : runtime_error(message){};

        ~Error() override = default;
    };

    virtual void run(std::string tped, std::string tfam,
                     std::vector<Position> &mutual_info,
                     size_t num_outputs) = 0;
};
#endif
