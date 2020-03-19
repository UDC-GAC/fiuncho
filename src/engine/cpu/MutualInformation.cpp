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
 * @file cpu/MutualInformation.cpp
 * @author Jorge González
 * @author Christian Ponte
 * @date 1 October 2018
 *
 * @brief MutualInformation class members implementation.
 */

#include "MutualInformation.h"
#include <algorithm>
#include <cfloat>
#include <cmath>

MutualInformation::MutualInformation(uint32_t numSNPs, uint16_t numCases,
                                     const std::vector<std::vector<uint32_t> *> &cases,
                                     uint16_t numCtrls, const std::vector<std::vector<uint32_t> *> &ctrls) :
        num_snps(numSNPs),
        cases(cases),
        num_cases(numCases),
        ctrls(ctrls),
        num_ctrls(numCtrls),
        _numEntriesCases(numCases / 32 + ((numCases % 32) > 0)),
        _numEntriesCtrls(numCtrls / 32 + ((numCtrls % 32) > 0)) {
    invInds = 1.0 / (numCases + numCtrls);

    float p = numCases * invInds;
    entY = (-1.0) * p * log2(p);

    p = numCtrls * invInds;
    entY -= p * log2(p);
}

long MutualInformation::compute(const std::vector<std::pair<uint32_t, uint32_t>> &pairs, uint16_t num_outputs,
                                Position *output) {
    uint32_t id1, id2, id3;
    float auxMIValue;
    long numAnal = 0;
    Position *auxMI;
    float minMI = FLT_MAX;
    uint16_t minMIPos = UINT16_MAX;
    uint16_t numEntriesWithMI = 0;
    DoubleContTable *doubleTable = new DoubleContTable(num_cases, num_ctrls);
    TripleContTable tripleTable;

    for (auto p : pairs) {
        //for (uint64_t iterPairs = 0; iterPairs < numPairs; iterPairs++) {

        id1 = p.first;
        id2 = p.second;

        _fillDoubleContTable(ctrls[id1], cases[id1], ctrls[id2], cases[id2], doubleTable);

        auxMIValue = _calcDoubleMI(doubleTable);

        // There are empty values in the array
        if (numEntriesWithMI < num_outputs) {
            auxMI = &output[numEntriesWithMI];
            auxMI->p1 = id1;
            auxMI->p2 = id2;
            auxMI->p3 = 0;
            auxMI->rank = auxMIValue;

            // If this is the minimum value of the array
            if (auxMIValue < minMI) {
                minMI = auxMIValue;
                minMIPos = numEntriesWithMI;
            }

            numEntriesWithMI++;
        } else if (auxMIValue > minMI) { // The value must be inserted
            auxMI = &output[minMIPos];
            auxMI->p1 = id1;
            auxMI->p2 = id2;
            auxMI->p3 = 0;
            auxMI->rank = auxMIValue;

            // Find the new minimum
            auxMI = std::min_element(output, output + num_outputs);
            minMI = auxMI->rank;
            uint16_t i = 0;
            while (1) {
                if (output[i].rank == minMI) {
                    break;
                }
                i++;
            }
            minMIPos = i;
        }

        for (id3 = id2 + 1; id3 < num_snps; id3++) {
            // Generate the contingency table of the 3-SNP
            _fillTripleContTable(doubleTable, &tripleTable, ctrls[id3], cases[id3]);

            // Calculate the MI with the contingency table
            auxMIValue = _calcTripleMI(&tripleTable);

            // There are empty values in the array
            if (numEntriesWithMI < num_outputs) {
                auxMI = &output[numEntriesWithMI];
                auxMI->p1 = id1;
                auxMI->p2 = id2;
                auxMI->p3 = id3;
                auxMI->rank = auxMIValue;

                // If this is the minimum value of the array
                if (auxMIValue < minMI) {
                    minMI = auxMIValue;
                    minMIPos = numEntriesWithMI;
                }

                numEntriesWithMI++;
            } else if (auxMIValue > minMI) { // The value must be inserted
                auxMI = &output[minMIPos];
                auxMI->p1 = id1;
                auxMI->p2 = id2;
                auxMI->p3 = id3;
                auxMI->rank = auxMIValue;

                // Find the new minimum
                auxMI = std::min_element(output, output + num_outputs);
                minMI = auxMI->rank;
                uint16_t i = 0;
                while (1) {
                    if (output[i].rank == minMI) {
                        break;
                    }
                    i++;
                }
                minMIPos = i;
            }

            numAnal++;
        }
    }

    delete doubleTable;
    return numAnal;
}

uint32_t popcount(uint32_t v) {
    uint32_t u;
    u = v - ((v >> 1) & 033333333333) - ((v >> 2) & 011111111111);
    return ((u + (u >> 3)) & 030707070707) % 63;
}

void MutualInformation::_fillDoubleContTable(std::vector<uint32_t> *s1_ctrls, std::vector<uint32_t> *s1_cases,
                                             std::vector<uint32_t> *s2_ctrls, std::vector<uint32_t> *s2_cases,
                                             DoubleContTable *table) {
    for (int i = 0; i < _numEntriesCases; i++) {
        table->_cases00[i] = s1_cases[0][i] & s2_cases[0][i];
        table->_cases01[i] = s1_cases[0][i] & s2_cases[1][i];
        table->_cases02[i] = s1_cases[0][i] & s2_cases[2][i];
        table->_cases10[i] = s1_cases[1][i] & s2_cases[0][i];
        table->_cases11[i] = s1_cases[1][i] & s2_cases[1][i];
        table->_cases12[i] = s1_cases[1][i] & s2_cases[2][i];
        table->_cases20[i] = s1_cases[2][i] & s2_cases[0][i];
        table->_cases21[i] = s1_cases[2][i] & s2_cases[1][i];
        table->_cases22[i] = s1_cases[2][i] & s2_cases[2][i];
    }

    for (int i = 0; i < _numEntriesCtrls; i++) {
        table->_ctrls00[i] = s1_ctrls[0][i] & s2_ctrls[0][i];
        table->_ctrls01[i] = s1_ctrls[0][i] & s2_ctrls[1][i];
        table->_ctrls02[i] = s1_ctrls[0][i] & s2_ctrls[2][i];
        table->_ctrls10[i] = s1_ctrls[1][i] & s2_ctrls[0][i];
        table->_ctrls11[i] = s1_ctrls[1][i] & s2_ctrls[1][i];
        table->_ctrls12[i] = s1_ctrls[1][i] & s2_ctrls[2][i];
        table->_ctrls20[i] = s1_ctrls[2][i] & s2_ctrls[0][i];
        table->_ctrls21[i] = s1_ctrls[2][i] & s2_ctrls[1][i];
        table->_ctrls22[i] = s1_ctrls[2][i] & s2_ctrls[2][i];
    }
}

float MutualInformation::_calcDoubleMI(DoubleContTable *table) {

    uint32_t contCases[9];
    uint32_t contCtrls[9];

    for (int i = 0; i < 9; i++) {
        contCases[i] = 0;
        contCtrls[i] = 0;
    }

    for (int i = 0; i < _numEntriesCases; i++) {
        contCases[0] += popcount(table->_cases00[i]);
        contCases[1] += popcount(table->_cases01[i]);
        contCases[2] += popcount(table->_cases02[i]);
        contCases[3] += popcount(table->_cases10[i]);
        contCases[4] += popcount(table->_cases11[i]);
        contCases[5] += popcount(table->_cases12[i]);
        contCases[6] += popcount(table->_cases20[i]);
        contCases[7] += popcount(table->_cases21[i]);
        contCases[8] += popcount(table->_cases22[i]);
    }
    for (int i = 0; i < _numEntriesCtrls; i++) {
        contCtrls[0] += popcount(table->_ctrls00[i]);
        contCtrls[1] += popcount(table->_ctrls01[i]);
        contCtrls[2] += popcount(table->_ctrls02[i]);
        contCtrls[3] += popcount(table->_ctrls10[i]);
        contCtrls[4] += popcount(table->_ctrls11[i]);
        contCtrls[5] += popcount(table->_ctrls12[i]);
        contCtrls[6] += popcount(table->_ctrls20[i]);
        contCtrls[7] += popcount(table->_ctrls21[i]);
        contCtrls[8] += popcount(table->_ctrls22[i]);
    }

    float entX = 0.0;
    float entAll = 0.0;
    float pCase, pCtrl;

    for (int i = 0; i < 9; i++) {
        pCase = contCases[i] * invInds;
        if (pCase != 0.0) {
            entAll -= pCase * log2(pCase);
        }

        pCtrl = contCtrls[i] * invInds;
        if (pCtrl != 0.0) {
            entAll -= pCtrl * log2(pCtrl);
        }

        pCase += pCtrl;
        if (pCase != 0.0) {
            entX -= pCase * log2(pCase);
        }
    }

    return entX + entY - entAll;
}

void MutualInformation::_fillTripleContTable(DoubleContTable *doubleTable, TripleContTable *tripleTable,
                                             std::vector<uint32_t> *s3_ctrls, std::vector<uint32_t> *s3_cases) {

    uint32_t aux;
    uint32_t auxSNP3Value;

    tripleTable->clearValues();

    for (int i = 0; i < _numEntriesCases; i++) {
        auxSNP3Value = s3_cases[0][i];

        aux = doubleTable->_cases00[i] & auxSNP3Value;
        tripleTable->_cases[0] += popcount(aux);

        aux = doubleTable->_cases01[i] & auxSNP3Value;
        tripleTable->_cases[1] += popcount(aux);

        aux = doubleTable->_cases02[i] & auxSNP3Value;
        tripleTable->_cases[2] += popcount(aux);

        aux = doubleTable->_cases10[i] & auxSNP3Value;
        tripleTable->_cases[3] += popcount(aux);

        aux = doubleTable->_cases11[i] & auxSNP3Value;
        tripleTable->_cases[4] += popcount(aux);

        aux = doubleTable->_cases12[i] & auxSNP3Value;
        tripleTable->_cases[5] += popcount(aux);

        aux = doubleTable->_cases20[i] & auxSNP3Value;
        tripleTable->_cases[6] += popcount(aux);

        aux = doubleTable->_cases21[i] & auxSNP3Value;
        tripleTable->_cases[7] += popcount(aux);

        aux = doubleTable->_cases22[i] & auxSNP3Value;
        tripleTable->_cases[8] += popcount(aux);


        auxSNP3Value = s3_cases[1][i];

        aux = doubleTable->_cases00[i] & auxSNP3Value;
        tripleTable->_cases[9] += popcount(aux);

        aux = doubleTable->_cases01[i] & auxSNP3Value;
        tripleTable->_cases[10] += popcount(aux);

        aux = doubleTable->_cases02[i] & auxSNP3Value;
        tripleTable->_cases[11] += popcount(aux);

        aux = doubleTable->_cases10[i] & auxSNP3Value;
        tripleTable->_cases[12] += popcount(aux);

        aux = doubleTable->_cases11[i] & auxSNP3Value;
        tripleTable->_cases[13] += popcount(aux);

        aux = doubleTable->_cases12[i] & auxSNP3Value;
        tripleTable->_cases[14] += popcount(aux);

        aux = doubleTable->_cases20[i] & auxSNP3Value;
        tripleTable->_cases[15] += popcount(aux);

        aux = doubleTable->_cases21[i] & auxSNP3Value;
        tripleTable->_cases[16] += popcount(aux);

        aux = doubleTable->_cases22[i] & auxSNP3Value;
        tripleTable->_cases[17] += popcount(aux);


        auxSNP3Value = s3_cases[2][i];

        aux = doubleTable->_cases00[i] & auxSNP3Value;
        tripleTable->_cases[18] += popcount(aux);

        aux = doubleTable->_cases01[i] & auxSNP3Value;
        tripleTable->_cases[19] += popcount(aux);

        aux = doubleTable->_cases02[i] & auxSNP3Value;
        tripleTable->_cases[20] += popcount(aux);

        aux = doubleTable->_cases10[i] & auxSNP3Value;
        tripleTable->_cases[21] += popcount(aux);

        aux = doubleTable->_cases11[i] & auxSNP3Value;
        tripleTable->_cases[22] += popcount(aux);

        aux = doubleTable->_cases12[i] & auxSNP3Value;
        tripleTable->_cases[23] += popcount(aux);

        aux = doubleTable->_cases20[i] & auxSNP3Value;
        tripleTable->_cases[24] += popcount(aux);

        aux = doubleTable->_cases21[i] & auxSNP3Value;
        tripleTable->_cases[25] += popcount(aux);

        aux = doubleTable->_cases22[i] & auxSNP3Value;
        tripleTable->_cases[26] += popcount(aux);
    }

    for (int i = 0; i < _numEntriesCtrls; i++) {
        auxSNP3Value = s3_ctrls[0][i];

        aux = doubleTable->_ctrls00[i] & auxSNP3Value;
        tripleTable->_ctrls[0] += popcount(aux);

        aux = doubleTable->_ctrls01[i] & auxSNP3Value;
        tripleTable->_ctrls[1] += popcount(aux);

        aux = doubleTable->_ctrls02[i] & auxSNP3Value;
        tripleTable->_ctrls[2] += popcount(aux);

        aux = doubleTable->_ctrls10[i] & auxSNP3Value;
        tripleTable->_ctrls[3] += popcount(aux);

        aux = doubleTable->_ctrls11[i] & auxSNP3Value;
        tripleTable->_ctrls[4] += popcount(aux);

        aux = doubleTable->_ctrls12[i] & auxSNP3Value;
        tripleTable->_ctrls[5] += popcount(aux);

        aux = doubleTable->_ctrls20[i] & auxSNP3Value;
        tripleTable->_ctrls[6] += popcount(aux);

        aux = doubleTable->_ctrls21[i] & auxSNP3Value;
        tripleTable->_ctrls[7] += popcount(aux);

        aux = doubleTable->_ctrls22[i] & auxSNP3Value;
        tripleTable->_ctrls[8] += popcount(aux);


        auxSNP3Value = s3_ctrls[1][i];

        aux = doubleTable->_ctrls00[i] & auxSNP3Value;
        tripleTable->_ctrls[9] += popcount(aux);

        aux = doubleTable->_ctrls01[i] & auxSNP3Value;
        tripleTable->_ctrls[10] += popcount(aux);

        aux = doubleTable->_ctrls02[i] & auxSNP3Value;
        tripleTable->_ctrls[11] += popcount(aux);

        aux = doubleTable->_ctrls10[i] & auxSNP3Value;
        tripleTable->_ctrls[12] += popcount(aux);

        aux = doubleTable->_ctrls11[i] & auxSNP3Value;
        tripleTable->_ctrls[13] += popcount(aux);

        aux = doubleTable->_ctrls12[i] & auxSNP3Value;
        tripleTable->_ctrls[14] += popcount(aux);

        aux = doubleTable->_ctrls20[i] & auxSNP3Value;
        tripleTable->_ctrls[15] += popcount(aux);

        aux = doubleTable->_ctrls21[i] & auxSNP3Value;
        tripleTable->_ctrls[16] += popcount(aux);

        aux = doubleTable->_ctrls22[i] & auxSNP3Value;
        tripleTable->_ctrls[17] += popcount(aux);


        auxSNP3Value = s3_ctrls[2][i];

        aux = doubleTable->_ctrls00[i] & auxSNP3Value;
        tripleTable->_ctrls[18] += popcount(aux);

        aux = doubleTable->_ctrls01[i] & auxSNP3Value;
        tripleTable->_ctrls[19] += popcount(aux);

        aux = doubleTable->_ctrls02[i] & auxSNP3Value;
        tripleTable->_ctrls[20] += popcount(aux);

        aux = doubleTable->_ctrls10[i] & auxSNP3Value;
        tripleTable->_ctrls[21] += popcount(aux);

        aux = doubleTable->_ctrls11[i] & auxSNP3Value;
        tripleTable->_ctrls[22] += popcount(aux);

        aux = doubleTable->_ctrls12[i] & auxSNP3Value;
        tripleTable->_ctrls[23] += popcount(aux);

        aux = doubleTable->_ctrls20[i] & auxSNP3Value;
        tripleTable->_ctrls[24] += popcount(aux);

        aux = doubleTable->_ctrls21[i] & auxSNP3Value;
        tripleTable->_ctrls[25] += popcount(aux);

        aux = doubleTable->_ctrls22[i] & auxSNP3Value;
        tripleTable->_ctrls[26] += popcount(aux);
    }
}

float MutualInformation::_calcTripleMI(TripleContTable *table) {

    float entX = 0.0;
    float entAll = 0.0;
    float pCase, pCtrl;

    for (int i = 0; i < 27; i++) {

        pCase = table->_cases[i] * invInds;
        if (pCase != 0.0) {
            entAll -= pCase * log2(pCase);
        }

        pCtrl = table->_ctrls[i] * invInds;
        if (pCtrl != 0.0) {
            entAll -= pCtrl * log2(pCtrl);
        }

        pCase += pCtrl;
        if (pCase != 0.0) {
            entX -= pCase * log2(pCase);
        }
    }

    return entX + entY - entAll;
}