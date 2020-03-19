/*
 * This file is part of MPI3SNP.
 * Copyright (C) 2014 - 2017 by Jorge González
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
 * @file gpu/EntropySearch.h
 * @author Jorge González
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief EntropySearch class declaration.
 */

#ifndef MPI3SNP_ENTROPYSEARCH_H
#define MPI3SNP_ENTROPYSEARCH_H

#include "GPUContTable.h"
#include <vector>
#include "Algorithm.h"

class MutualInformation : public Algorithm<uint2> {
public:
    ~MutualInformation();

    MutualInformation(bool isMI, uint32_t numSNPs, uint16_t numCases, uint16_t numCtrls,
                  std::vector<std::vector<uint32_t> *> cases, std::vector<std::vector<uint32_t> *> ctrls);

    long compute(const std::vector<uint2> &workload, uint16_t num_outputs, Position *output) override;

private:
    bool _isMI;

    uint32_t _numSNPs;
    uint16_t _numEntriesCase;
    uint16_t _numEntriesCtrl;

    uint32_t *_dev0Cases;
    uint32_t *_dev1Cases;
    uint32_t *_dev2Cases;
    uint32_t *_dev0Ctrls;
    uint32_t *_dev1Ctrls;
    uint32_t *_dev2Ctrls;
    uint2 *_devIds;

    void _findNHighestMI(uint3 *_hostMiIds, float *_hostMIValues, uint64_t totalValues, float &minMI, uint16_t &minMIPos, uint16_t &numEntriesWithMI,
                         size_t num_outputs, Position *output);
};


#endif //MPI3SNP_ENTROPYSEARCH_H
