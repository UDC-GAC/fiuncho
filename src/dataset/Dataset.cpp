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
 * @brief Implementation of Dataset class.
 */

#include "Dataset.h"
#include <fstream>
#include <functional>
#include <numeric>

Dataset::Dataset(std::string tped_path, std::string tfam_path, Representation rep) {
    std::ifstream file;
    std::vector<Individual> individuals;
    std::vector<SNP> snps;

    file.open(tfam_path.c_str(), std::ios::in);
    if (!file.is_open()) {
        throw Read_error("Error while opening " + tfam_path + ", check file path/permissions");
    }
    try {
        Individual ind;
        while (file >> ind) {
            individuals.push_back(ind);
        }
    } catch (const Individual::InvalidIndividual &e) {
        throw Read_error("Error in " + tfam_path + ":" + std::to_string(individuals.size() + 1) + ": " + e.what());
    }
    file.close();

    file.open(tped_path.c_str(), std::ios::in);
    if (!file.is_open()) {
        throw Read_error("Error while opening " + tped_path + ", check file path/permissions");
    }
    try {
        SNP snp;
        while (file >> snp) {
            if (snp.genotypes.size() == individuals.size()) {
                snps.push_back(snp);
            } else {
                throw Read_error("Error in " + tped_path + ":" + std::to_string(snps.size() + 1) +
                                ": the number of nucleotides does not match the number of individuals");
            }
        }
    } catch (const SNP::InvalidSNP &e) {
        throw Read_error("Error in " + tped_path + ":" + std::to_string(snps.size() + 1) + ": " + e.what());
    }
    file.close();

    snp_count = snps.size();

    switch (rep) {
        case Regular:
            regular_representation(individuals, snps);
            break;
        case Transposed:
            transposed_representation(individuals, snps);
            break;
    }
}

unsigned long find_index(unsigned long start, unsigned long end, std::function<bool(unsigned long)> fun) {
    while (start != end && !fun(start)) {
        start++;
    }
    return start;
}

void Dataset::regular_representation(std::vector<Individual> &inds, std::vector<SNP> &snps) {
    std::vector<uint32_t> *cases_vec = nullptr;
    std::vector<uint32_t> *ctrls_vec = nullptr;
    uint32_t cases_buff[3], ctrls_buff[3];
    unsigned int n_cases_buff = 0, n_ctrls_buff = 0;
    size_t i, j;

    for (SNP s : snps) {
        // Initialize buffers
        n_ctrls_buff = 0;
        n_cases_buff = 0;
        for (j = 0; j < 3; j++) {
            ctrls_buff[j] = 0;
            cases_buff[j] = 0;
        }
        ctrls_vec = new std::vector<uint32_t>[3];
        cases_vec = new std::vector<uint32_t>[3];
        for (i = 0; i < inds.size(); i++) {
            if (inds[i].ph == 1) {
                for (j = 0; j < 3; j++) {
                    ctrls_buff[j] <<= 1;
                    ctrls_buff[j] += s.genotypes[i] == j;
                }
                if (++n_ctrls_buff == 32) {
                    n_ctrls_buff = 0;
                    for (j = 0; j < 3; j++) {
                        ctrls_vec[j].push_back(ctrls_buff[j]);
                        ctrls_buff[j] = 0;
                    }
                }
            } else {
                for (j = 0; j < 3; j++) {
                    cases_buff[j] <<= 1;
                    cases_buff[j] += s.genotypes[i] == j;
                }
                if (++n_cases_buff == 32) {
                    n_cases_buff = 0;
                    for (j = 0; j < 3; j++) {
                        cases_vec[j].push_back(cases_buff[j]);
                        cases_buff[j] = 0;
                    }
                }
            }
        }
        if (n_ctrls_buff > 0) {
            for (j = 0; j < 3; j++) {
                ctrls_vec[j].push_back(ctrls_buff[j]);
            }
        }
        if (n_cases_buff > 0) {
            for (j = 0; j < 3; j++) {
                cases_vec[j].push_back(cases_buff[j]);
            }
        }
        ctrls.push_back(ctrls_vec);
        cases.push_back(cases_vec);
    }

    num_ctrls = (ctrls[0][0].size() - 1) * 32 + n_ctrls_buff;
    num_cases = (cases[0][0].size() - 1) * 32 + n_cases_buff;
}

void Dataset::transposed_representation(std::vector<Individual> &inds, std::vector<SNP> &snps) {
    std::vector<unsigned long> scases(32), sctrls(32);
    unsigned long ctrlpos = 0, casepos = 0;
    uint32_t cases_buffer[3], ctrls_buffer[3];
    size_t i, j;

    std::vector<uint32_t> *t_cases = new std::vector<uint32_t>[3];
    std::vector<uint32_t> *t_ctrls = new std::vector<uint32_t>[3];

    num_cases = 0;
    num_ctrls = 0;

    // Iterate on all SNPs considering a block of 32 cases and controls, for all individuals
    while (num_cases + num_ctrls < inds.size()) {
        // Select next 32 controls and cases
        scases.clear();
        while (scases.size() < 32 &&
               (casepos = find_index(casepos, inds.size(), [&inds](int k) { return inds[k].ph == 2; })) < inds.size()) {
            scases.push_back(casepos++);
        }
        sctrls.clear();
        while (sctrls.size() < 32 &&
               (ctrlpos = find_index(ctrlpos, inds.size(), [&inds](int k) { return inds[k].ph == 1; })) < inds.size()) {
            sctrls.push_back(ctrlpos++);
        }
        // Read all SNPs for those 32 controls and cases
        for (i = 0; i < snps.size(); i++) {
            for (j = 0; j < 3; j++) {
                cases_buffer[j] = 0;
                ctrls_buffer[j] = 0;
            }
            for (unsigned long pos : scases) {
                for (j = 0; j < 3; j++) {
                    cases_buffer[j] = cases_buffer[j] << 1;
                    cases_buffer[j] += snps[i].genotypes[pos] == j;
                }
            }
            for (unsigned long pos : sctrls) {
                for (j = 0; j < 3; j++) {
                    ctrls_buffer[j] = ctrls_buffer[j] << 1;
                    ctrls_buffer[j] += snps[i].genotypes[pos] == j;
                }
            }
            // Save buffers when not empty
            if (!scases.empty()) {
                for (j = 0; j < 3; j++) {
                    t_cases[j].push_back(cases_buffer[j]);
                }
            }
            if (!sctrls.empty()) {
                for (j = 0; j < 3; j++) {
                    t_ctrls[j].push_back(ctrls_buffer[j]);
                }
            }
        }
        num_cases += scases.size();
        num_ctrls += sctrls.size();
    }

    this->cases.push_back(t_cases);
    this->ctrls.push_back(t_ctrls);
}

Dataset::~Dataset() {
    for (std::vector<uint32_t> *item : cases) {
        delete[] item;
    }
    for (std::vector<uint32_t> *item : ctrls) {
        delete[] item;
    }
}