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
 * @file MutualInformation.h
 * @author Christian Ponte
 * @date 30 May 2018
 *
 * @brief MutualInformation class declaration.
 */

#ifndef FIUNCHO_MUTUALINFORMATION_H
#define FIUNCHO_MUTUALINFORMATION_H

#include <cstdint>
#include <fiuncho/engine/Algorithm.h>
#include <fiuncho/engine/cpu/ContTable.h>
#include <vector>

class MutualInformation : public Algorithm<std::pair<uint32_t, uint32_t>> {
  public:
    MutualInformation(uint32_t numSNPs, uint16_t numCases,
                      const std::vector<std::vector<uint32_t> *> &cases,
                      uint16_t numCtrls,
                      const std::vector<std::vector<uint32_t> *> &ctrls);

    long compute(const std::vector<std::pair<uint32_t, uint32_t>> &pairs,
                 uint16_t num_outputs, Position *output) override;

  private:
    void _fillDoubleContTable(std::vector<uint32_t> *s1_ctrls,
                              std::vector<uint32_t> *s1_cases,
                              std::vector<uint32_t> *s2_ctrls,
                              std::vector<uint32_t> *s2_cases,
                              DoubleContTable *table);

    float _calcDoubleMI(DoubleContTable *table);

    void _fillTripleContTable(DoubleContTable *doubleTable,
                              TripleContTable *tripleTable,
                              std::vector<uint32_t> *s3_ctrls,
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

#endif
