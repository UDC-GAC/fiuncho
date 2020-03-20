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
 * @file Dataset.h
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief Distributor class declaration. This class implements the parallel
 * distribution strategy, shared by both cpu and gpu implementation.
 */

#ifndef FIUNCHO_DISTRIBUTOR_H
#define FIUNCHO_DISTRIBUTOR_H

#include <vector>

template <typename T, typename U> class Distributor {
  public:
    Distributor(const T &size, const unsigned int &frac)
        : size(size), frac(frac){};

    void get_pairs(U (*constructor)(T, T), const unsigned int &id,
                   std::vector<U> &values) {
        const size_t total_pairs = size * (size - 1) / 2;
        values.reserve(total_pairs / frac + (id < total_pairs % frac));
        T i, j, offset;
        offset = id;
        for (i = 0; i < size; i++) {
            for (j = i + 1 + offset; j < size; j += frac) {
                values.push_back(constructor(i, j));
            }
            offset = j - size;
        }
    }

  private:
    const T size;
    const unsigned int frac;
};

#endif
