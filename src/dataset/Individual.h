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
 * @file Individual.h
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief Individual class definition and implementation.
 */

#ifndef MPI3SNP_INDIVIDUAL_H
#define MPI3SNP_INDIVIDUAL_H

#include <string>
#include <sstream>

struct Individual {
    class InvalidIndividual : public std::runtime_error {
    public:
        InvalidIndividual(const std::string &message) : std::runtime_error(message) {};

        virtual ~InvalidIndividual() {};
    };

    std::string fid; // Family ID
    std::string iid; // Within-family ID (cannot be '0')
    std::string f_iid; // Within-family ID of father ('0' if father isn't in dataset)
    std::string m_iid; // Within-family ID of mother ('0' if mother isn't in dataset)
    int sex; // Sex code ('1' = male, '2' = female, '0' = unknown)
    int ph; // Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)

    friend std::istream &operator>>(std::istream &str, Individual &ind) {
        std::string line;
        Individual tmp;
        if (std::getline(str, line)) {
            std::stringstream iss(line);
            if (iss >> tmp.fid &&
                iss >> tmp.iid &&
                iss >> tmp.f_iid &&
                iss >> tmp.m_iid &&
                iss >> tmp.sex &&
                iss >> tmp.ph && (tmp.ph == 1 || tmp.ph == 2)) {
                /* OK: All read operations worked */
                ind.swap(tmp);  // C++03 as this answer was written a long time ago.
            } else {
                throw InvalidIndividual("parsing error at " + std::to_string(iss.tellg()));
            }
        }
        return str;
    }

    void swap(Individual &other) {
        std::swap(fid, other.fid);
        std::swap(iid, other.iid);
        std::swap(f_iid, other.f_iid);
        std::swap(m_iid, other.m_iid);
        std::swap(sex, other.sex);
        std::swap(ph, other.ph);
    }
};

#endif //MPI3SNP_INDIVIDUAL_H
