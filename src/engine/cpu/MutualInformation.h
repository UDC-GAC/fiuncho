/*
 * This file is part of MPI3SNP.
 * Copyright (C) 2018 by Christian Ponte
 *
 * MPI3SNP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MPI3SNP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MPI3SNP. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file cpu/MutualInformation.h
 * @author Jorge González
 * @author Christian Ponte
 * @date 1 October 2018
 *
 * @brief MutualInformation class declaration.
 */

#ifndef MPI3SNP_MUTUALINFORMATION_H
#define MPI3SNP_MUTUALINFORMATION_H

#include <cstdint>
#include <vector>
#include "Algorithm.h"
#include "ContTable.h"

class MutualInformation : public Algorithm<std::pair<uint32_t, uint32_t>> {
public:
    MutualInformation(uint32_t numSNPs, uint16_t numCases, const std::vector<std::vector<uint32_t> *> &cases,
                      uint16_t numCtrls, const std::vector<std::vector<uint32_t> *> &ctrls);

    long compute(const std::vector<std::pair<uint32_t, uint32_t>> &pairs, uint16_t num_outputs, Position *output)
    override;

private:

    void _fillDoubleContTable(std::vector<uint32_t> *s1_ctrls, std::vector<uint32_t> *s1_cases,
                              std::vector<uint32_t> *s2_ctrls, std::vector<uint32_t> *s2_cases, DoubleContTable *table);

    float _calcDoubleMI(DoubleContTable *table);

    void
    _fillTripleContTable(DoubleContTable *doubleTable, TripleContTable *tripleTable, std::vector<uint32_t> *s3_ctrls,
                         std::vector<uint32_t> *s3_cases);

    float _calcTripleMI(TripleContTable *table);

    const uint32_t num_snps;
    const std::vector<std::vector<uint32_t> *> &cases;
    const uint16_t num_cases;
    const std::vector<std::vector<uint32_t> *> &ctrls;
    const uint16_t num_ctrls;
    const uint16_t _numEntriesCases;
    const uint16_t _numEntriesCtrls;
    float invInds;
    // Entropy of Y
    float entY;
};


#endif //MPI3SNP_ENTROPYSEARCH_H
