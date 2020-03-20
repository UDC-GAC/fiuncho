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

#ifndef FIUNCHO_POSITION_H
#define FIUNCHO_POSITION_H

#include <cinttypes>
#include <string>

class Position {
  public:
    uint32_t p1, p2, p3;
    float rank;

    bool operator<(const Position p) const { return rank < p.rank; }

    inline std::string to_string() {
        return std::to_string(p1) + " " + std::to_string(p2) + " " +
               std::to_string(p3) + " " + std::to_string(rank);
    }
};

#endif
