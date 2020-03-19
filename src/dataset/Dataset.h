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
 * @brief Dataset class declaration, responsible for reading the data and storing it using a bitwise representation.
 */

#ifndef MPI3SNP_DATASET_H
#define MPI3SNP_DATASET_H

#include <string>
#include <stdexcept>
#include <vector>
#include "Individual.h"
#include "SNP.h"

class Dataset {
public:
    class Read_error : public std::runtime_error {
    public:
        Read_error(const std::string message) : runtime_error(message) {}

        virtual ~Read_error() {};
    };

    enum Representation {
        Regular,
        Transposed
    };

    Dataset(std::string tped_path, std::string tfam_path, Representation rep);

    ~Dataset();

    inline std::vector<std::vector<uint32_t> *> &get_cases() {
        return cases;
    }

    inline std::vector<std::vector<uint32_t> *> &get_ctrls() {
        return ctrls;
    }

    inline uint32_t get_SNP_count() {
        return snp_count;
    }

    inline uint16_t get_ctrl_count() {
        return num_ctrls;
    }

    inline uint16_t get_case_count() {
        return num_cases;
    }

private:
    void regular_representation(std::vector<Individual> &inds, std::vector<SNP> &snps);

    void transposed_representation(std::vector<Individual> &inds, std::vector<SNP> &snps);

    std::vector<std::vector<uint32_t> *> cases;
    std::vector<std::vector<uint32_t> *> ctrls;
    uint32_t snp_count;
    uint16_t num_cases, num_ctrls;
};

#endif //MPI3SNP_DATASET_H
