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
#include <memory>

template <typename T> class MaxArray
{
    const size_t maxsize;
    std::unique_ptr<T[]> ptr;
    T *a;
    // The position of the minimum value
    T *min_pos;
    // Number of entries of the array full
    size_t current_size;

  public:
    MaxArray(const MaxArray<T> &) = delete;
    MaxArray(MaxArray<T> &&) = default;

    explicit MaxArray(const size_t &maxsize)
        : maxsize(maxsize), ptr(new T[maxsize]), a(ptr.get()), min_pos(a),
          current_size(0)
    {
    }

    T &operator[](std::size_t pos) { return a[pos]; }

    void add(const T &value)
    {
        // There are empty values in the array
        if (current_size < maxsize) {
            T *pos = &a[current_size++];
            *pos = value;
            // If this is the minimum value of the array
            if (*pos < *min_pos) {
                min_pos = pos;
            }
        } else if (value > *min_pos) { // The value must be inserted
            *min_pos = value;
            // Find the new minimum
            min_pos = std::min_element(a, a + maxsize);
        }
    }

    T *array()
    {
        auto copy = new T[maxsize];
        memcpy(copy, a, sizeof(T) * maxsize);
        return copy;
    }

    size_t size() const { return current_size; }
};

#endif
