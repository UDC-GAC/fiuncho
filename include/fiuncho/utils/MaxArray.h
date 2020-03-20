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

#ifndef FIUNCHO_MAXARRAY_H
#define FIUNCHO_MAXARRAY_H

#include <algorithm>
#include <cstdlib>
#include <cstring>

template <typename T> class MaxArray {
  public:
    explicit MaxArray(const size_t &size)
        : size(size), a(new T[size]), min_pos(a), emplacement(0) {}

    ~MaxArray() { delete[](a); }

    T &operator[](std::size_t pos) { return a[pos]; }

    void add(const T &value) {
        // There are empty values in the array
        if (emplacement < size) {
            T *pos = &a[emplacement++];
            *pos = value;
            // If this is the minimum value of the array
            if (*pos < *min_pos) {
                min_pos = pos;
            }
        } else if (value > *min_pos) { // The value must be inserted
            *min_pos = value;
            // Find the new minimum
            min_pos = std::min_element(a, a + size);
        }
    }

    T *array() {
        auto copy = new T[size];
        memcpy(copy, a, sizeof(T) * size);
        return copy;
    }

  private:
    T *a;
    // The position of the minimum value
    T *min_pos;
    // Number of entries of the array full
    size_t emplacement;

    const size_t size;
};

#endif
