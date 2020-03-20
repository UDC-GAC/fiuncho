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
 * @date 29 May 2018
 *
 * @brief Dataset class declaration, responsible for reading the data and
 * storing it using a bitwise representation.
 */

#ifndef FIUNCHO_DATASET_H
#define FIUNCHO_DATASET_H

#include <fiuncho/dataset/Individual.h>
#include <fiuncho/dataset/SNP.h>
#include <string>
#include <vector>

class Dataset {
  public:
    class Read_error : public std::runtime_error {
      public:
        Read_error(const std::string message) : runtime_error(message) {}

        virtual ~Read_error(){};
    };

    enum Representation { Regular, Transposed };

    Dataset(std::string tped_path, std::string tfam_path, Representation rep);

    ~Dataset();

    inline std::vector<std::vector<uint32_t> *> &get_cases() { return cases; }

    inline std::vector<std::vector<uint32_t> *> &get_ctrls() { return ctrls; }

    inline uint32_t get_SNP_count() { return snp_count; }

    inline uint16_t get_ctrl_count() { return num_ctrls; }

    inline uint16_t get_case_count() { return num_cases; }

  private:
    void regular_representation(std::vector<Individual> &inds,
                                std::vector<SNP> &snps);

    void transposed_representation(std::vector<Individual> &inds,
                                   std::vector<SNP> &snps);

    std::vector<std::vector<uint32_t> *> cases;
    std::vector<std::vector<uint32_t> *> ctrls;
    uint32_t snp_count;
    uint16_t num_cases, num_ctrls;
};

#endif
