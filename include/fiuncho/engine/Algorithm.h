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

#ifndef FIUNCHO_ALGORITHM_H
#define FIUNCHO_ALGORITHM_H

#include <fiuncho/engine/Position.h>
#include <vector>

template <typename T> class Algorithm {
  public:
    virtual ~Algorithm() = default;

    virtual long compute(const std::vector<T> &workload, uint16_t num_outputs,
                         Position *output) = 0;
};

#endif
