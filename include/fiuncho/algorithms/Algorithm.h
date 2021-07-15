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

#include <fiuncho/ContingencyTable.h>

template <class T>
class Algorithm {
  public:
    virtual ~Algorithm() = default;

    template <class U>
    T compute(const ContingencyTable<U> &table) const noexcept;
};

#endif
