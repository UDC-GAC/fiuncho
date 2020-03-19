/*
 * MPI3SNP ~ https://github.com/chponte/mpi3snp
 *
 * Copyright 2018 Christian Ponte
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
 * documentation files (the "Software"), to deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the
 * Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
 * WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * @file Dataset.h
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief Distributor class declaration. This class implements the parallel distribution strategy, shared by both cpu
 * and gpu implementation.
 */

#ifndef MPI3SNP_DISTRIBUTOR_H
#define MPI3SNP_DISTRIBUTOR_H

#include <vector>

template<typename T, typename U>
class Distributor {
public:
    Distributor(const T &size, const unsigned int &frac) : size(size), frac(frac) {};

    void get_pairs(U (*constructor)(T, T), const unsigned int &id, std::vector<U> &values) {
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

#endif //MPI3SNP_DISTRIBUTOR_H
