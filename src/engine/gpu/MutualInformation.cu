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
 * @file gpu/EntropySearch.cu
 * @author Jorge González
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief EntropySearch class members implementation.
 */

#include "MutualInformation.h"
#include <float.h>
#include <algorithm>

/*number of streaming processor scale*/
#define NUM_TH_PER_BLOCK    16

__constant__ float _invInds;
__constant__ float _entY;
__constant__ float _MAX_FLOAT;

static __device__ uint32_t
_devpopcount(uint32_t
v) {
//uint32_t u;
//u = v - ((v>>1) & 033333333333) - ((v>>2) & 011111111111);
// return ((u + (u>>3)) & 030707070707) % 63;
return
__popc(v);
}

static __global__ void _kernelDoubleTable(uint64_t numPairs, uint32_t numSNPs, uint16_t numEntriesCases,
                                          uint16_t numEntriesCtrls, uint2 *devIds, uint32_t *dev0Cases,
                                          uint32_t *dev1Cases, uint32_t *dev2Cases, uint32_t *dev0Ctrls,
                                          uint32_t *dev1Ctrls, uint32_t *dev2Ctrls,
                                          GPUDoubleContTable *doubleTables) {
    int gid = threadIdx.x + blockIdx.x * blockDim.x; /*global id*/

    if (gid >= numPairs) {
        return;
    }

    uint32_t
    myId1 = devIds[gid].x;
    uint32_t
    myId2 = devIds[gid].y;
    uint32_t
    entry1 = myId1;
    uint32_t
    entry2 = myId2;

    GPUDoubleContTable *table = &doubleTables[gid];

    for (int i = 0; i < numEntriesCases; i++, entry1 += numSNPs, entry2 += numSNPs) {
        table->_cases00[i] = dev0Cases[entry1] & dev0Cases[entry2];
        table->_cases01[i] = dev0Cases[entry1] & dev1Cases[entry2];
        table->_cases02[i] = dev0Cases[entry1] & dev2Cases[entry2];
        table->_cases10[i] = dev1Cases[entry1] & dev0Cases[entry2];
        table->_cases11[i] = dev1Cases[entry1] & dev1Cases[entry2];
        table->_cases12[i] = dev1Cases[entry1] & dev2Cases[entry2];
        table->_cases20[i] = dev2Cases[entry1] & dev0Cases[entry2];
        table->_cases21[i] = dev2Cases[entry1] & dev1Cases[entry2];
        table->_cases22[i] = dev2Cases[entry1] & dev2Cases[entry2];
    }

    entry1 = myId1;
    entry2 = myId2;

    for (int i = 0; i < numEntriesCtrls; i++, entry1 += numSNPs, entry2 += numSNPs) {
        table->_ctrls00[i] = dev0Ctrls[entry1] & dev0Ctrls[entry2];
        table->_ctrls01[i] = dev0Ctrls[entry1] & dev1Ctrls[entry2];
        table->_ctrls02[i] = dev0Ctrls[entry1] & dev2Ctrls[entry2];
        table->_ctrls10[i] = dev1Ctrls[entry1] & dev0Ctrls[entry2];
        table->_ctrls11[i] = dev1Ctrls[entry1] & dev1Ctrls[entry2];
        table->_ctrls12[i] = dev1Ctrls[entry1] & dev2Ctrls[entry2];
        table->_ctrls20[i] = dev2Ctrls[entry1] & dev0Ctrls[entry2];
        table->_ctrls21[i] = dev2Ctrls[entry1] & dev1Ctrls[entry2];
        table->_ctrls22[i] = dev2Ctrls[entry1] & dev2Ctrls[entry2];
    }
}

static __global__ void _kernelTripleMI(uint64_t numPairs, uint32_t numSNPs,
                                       uint16_t numEntriesCases, uint16_t numEntriesCtrls, uint2 *devIds,
                                       uint32_t *dev0Cases, uint32_t *dev1Cases, uint32_t *dev2Cases,
                                       uint32_t *dev0Ctrls, uint32_t *dev1Ctrls, uint32_t *dev2Ctrls,
                                       GPUDoubleContTable *devDoubleTables, uint16_t numOutputs,
                                       float *devMIValues, uint3 *devMiIds) {

    extern __shared__ uint32_t
    sharedMem[];
    uint32_t * cases00 = sharedMem;
    uint32_t * cases01 = &sharedMem[numEntriesCases];
    uint32_t * cases02 = &sharedMem[2 * numEntriesCases];
    uint32_t * cases10 = &sharedMem[3 * numEntriesCases];
    uint32_t * cases11 = &sharedMem[4 * numEntriesCases];
    uint32_t * cases12 = &sharedMem[5 * numEntriesCases];
    uint32_t * cases20 = &sharedMem[6 * numEntriesCases];
    uint32_t * cases21 = &sharedMem[7 * numEntriesCases];
    uint32_t * cases22 = &sharedMem[8 * numEntriesCases];
    uint32_t * ctrls00 = &sharedMem[9 * numEntriesCases];
    uint32_t * ctrls01 = &sharedMem[9 * numEntriesCases + numEntriesCtrls];
    uint32_t * ctrls02 = &sharedMem[9 * numEntriesCases + 2 * numEntriesCtrls];
    uint32_t * ctrls10 = &sharedMem[9 * numEntriesCases + 3 * numEntriesCtrls];
    uint32_t * ctrls11 = &sharedMem[9 * numEntriesCases + 4 * numEntriesCtrls];
    uint32_t * ctrls12 = &sharedMem[9 * numEntriesCases + 5 * numEntriesCtrls];
    uint32_t * ctrls20 = &sharedMem[9 * numEntriesCases + 6 * numEntriesCtrls];
    uint32_t * ctrls21 = &sharedMem[9 * numEntriesCases + 7 * numEntriesCtrls];
    uint32_t * ctrls22 = &sharedMem[9 * numEntriesCases + 8 * numEntriesCtrls];
    uint32_t * shMIId = &sharedMem[9 * (numEntriesCases + numEntriesCtrls)];
    float *shMIValues = (float *) &sharedMem[9 * (numEntriesCases + numEntriesCtrls) + blockDim.x * numOutputs];

    uint32_t * myOutIds = &shMIId[threadIdx.x * numOutputs];
    float *myOutValues = &shMIValues[threadIdx.x * numOutputs];

    uint16_t numEntriesWithMI = 0;
    float minMI = _MAX_FLOAT;
    uint16_t minMIPos = 0;

    // Copy the information of the double contingency table to shared memory
    GPUDoubleContTable *table = &devDoubleTables[blockIdx.x];

    for (int iter = threadIdx.x; iter < numEntriesCases; iter += blockDim.x) {
        cases00[iter] = table->_cases00[iter];
        cases01[iter] = table->_cases01[iter];
        cases02[iter] = table->_cases02[iter];
        cases10[iter] = table->_cases10[iter];
        cases11[iter] = table->_cases11[iter];
        cases12[iter] = table->_cases12[iter];
        cases20[iter] = table->_cases20[iter];
        cases21[iter] = table->_cases21[iter];
        cases22[iter] = table->_cases22[iter];
    }
    for (int iter = threadIdx.x; iter < numEntriesCtrls; iter += blockDim.x) {
        ctrls00[iter] = table->_ctrls00[iter];
        ctrls01[iter] = table->_ctrls01[iter];
        ctrls02[iter] = table->_ctrls02[iter];
        ctrls10[iter] = table->_ctrls10[iter];
        ctrls11[iter] = table->_ctrls11[iter];
        ctrls12[iter] = table->_ctrls12[iter];
        ctrls20[iter] = table->_ctrls20[iter];
        ctrls21[iter] = table->_ctrls21[iter];
        ctrls22[iter] = table->_ctrls22[iter];
    }

    __syncthreads();

    uint16_t tripleCases[27];
    uint16_t tripleCtrls[27];

    uint32_t
    myId1 = devIds[blockIdx.x].x;
    uint32_t
    myId2 = devIds[blockIdx.x].y;
    uint32_t
    iterId3 = myId2 + threadIdx.x + 1;

    uint32_t
    aux;
    uint32_t
    auxSNP3Value;

    // Each thread computes several triples using the same pair
    for (; iterId3 < numSNPs; iterId3 += blockDim.x) {
        // Starts creating the contingency table
        for (int i = 0; i < 27; i++) {
            tripleCases[i] = 0;
            tripleCtrls[i] = 0;
        }

        for (int i = 0; i < numEntriesCases; i++) {
            auxSNP3Value = dev0Cases[iterId3 + i * numSNPs];

            aux = cases00[i] & auxSNP3Value;
            tripleCases[0] += _devpopcount(aux);

            aux = cases01[i] & auxSNP3Value;
            tripleCases[1] += _devpopcount(aux);

            aux = cases02[i] & auxSNP3Value;
            tripleCases[2] += _devpopcount(aux);

            aux = cases10[i] & auxSNP3Value;
            tripleCases[3] += _devpopcount(aux);

            aux = cases11[i] & auxSNP3Value;
            tripleCases[4] += _devpopcount(aux);

            aux = cases12[i] & auxSNP3Value;
            tripleCases[5] += _devpopcount(aux);

            aux = cases20[i] & auxSNP3Value;
            tripleCases[6] += _devpopcount(aux);

            aux = cases21[i] & auxSNP3Value;
            tripleCases[7] += _devpopcount(aux);

            aux = cases22[i] & auxSNP3Value;
            tripleCases[8] += _devpopcount(aux);

            auxSNP3Value = dev1Cases[iterId3 + i * numSNPs];

            aux = cases00[i] & auxSNP3Value;
            tripleCases[9] += _devpopcount(aux);

            aux = cases01[i] & auxSNP3Value;
            tripleCases[10] += _devpopcount(aux);

            aux = cases02[i] & auxSNP3Value;
            tripleCases[11] += _devpopcount(aux);

            aux = cases10[i] & auxSNP3Value;
            tripleCases[12] += _devpopcount(aux);

            aux = cases11[i] & auxSNP3Value;
            tripleCases[13] += _devpopcount(aux);

            aux = cases12[i] & auxSNP3Value;
            tripleCases[14] += _devpopcount(aux);

            aux = cases20[i] & auxSNP3Value;
            tripleCases[15] += _devpopcount(aux);

            aux = cases21[i] & auxSNP3Value;
            tripleCases[16] += _devpopcount(aux);

            aux = cases22[i] & auxSNP3Value;
            tripleCases[17] += _devpopcount(aux);

            auxSNP3Value = dev2Cases[iterId3 + i * numSNPs];

            aux = cases00[i] & auxSNP3Value;
            tripleCases[18] += _devpopcount(aux);

            aux = cases01[i] & auxSNP3Value;
            tripleCases[19] += _devpopcount(aux);

            aux = cases02[i] & auxSNP3Value;
            tripleCases[20] += _devpopcount(aux);

            aux = cases10[i] & auxSNP3Value;
            tripleCases[21] += _devpopcount(aux);

            aux = cases11[i] & auxSNP3Value;
            tripleCases[22] += _devpopcount(aux);

            aux = cases12[i] & auxSNP3Value;
            tripleCases[23] += _devpopcount(aux);

            aux = cases20[i] & auxSNP3Value;
            tripleCases[24] += _devpopcount(aux);

            aux = cases21[i] & auxSNP3Value;
            tripleCases[25] += _devpopcount(aux);

            aux = cases22[i] & auxSNP3Value;
            tripleCases[26] += _devpopcount(aux);
        }

        for (int i = 0; i < numEntriesCtrls; i++) {
            auxSNP3Value = dev0Ctrls[iterId3 + i * numSNPs];

            aux = ctrls00[i] & auxSNP3Value;
            tripleCtrls[0] += _devpopcount(aux);

            aux = ctrls01[i] & auxSNP3Value;
            tripleCtrls[1] += _devpopcount(aux);

            aux = ctrls02[i] & auxSNP3Value;
            tripleCtrls[2] += _devpopcount(aux);

            aux = ctrls10[i] & auxSNP3Value;
            tripleCtrls[3] += _devpopcount(aux);

            aux = ctrls11[i] & auxSNP3Value;
            tripleCtrls[4] += _devpopcount(aux);

            aux = ctrls12[i] & auxSNP3Value;
            tripleCtrls[5] += _devpopcount(aux);

            aux = ctrls20[i] & auxSNP3Value;
            tripleCtrls[6] += _devpopcount(aux);

            aux = ctrls21[i] & auxSNP3Value;
            tripleCtrls[7] += _devpopcount(aux);

            aux = ctrls22[i] & auxSNP3Value;
            tripleCtrls[8] += _devpopcount(aux);

            auxSNP3Value = dev1Ctrls[iterId3 + i * numSNPs];

            aux = ctrls00[i] & auxSNP3Value;
            tripleCtrls[9] += _devpopcount(aux);

            aux = ctrls01[i] & auxSNP3Value;
            tripleCtrls[10] += _devpopcount(aux);

            aux = ctrls02[i] & auxSNP3Value;
            tripleCtrls[11] += _devpopcount(aux);

            aux = ctrls10[i] & auxSNP3Value;
            tripleCtrls[12] += _devpopcount(aux);

            aux = ctrls11[i] & auxSNP3Value;
            tripleCtrls[13] += _devpopcount(aux);

            aux = ctrls12[i] & auxSNP3Value;
            tripleCtrls[14] += _devpopcount(aux);

            aux = ctrls20[i] & auxSNP3Value;
            tripleCtrls[15] += _devpopcount(aux);

            aux = ctrls21[i] & auxSNP3Value;
            tripleCtrls[16] += _devpopcount(aux);

            aux = ctrls22[i] & auxSNP3Value;
            tripleCtrls[17] += _devpopcount(aux);

            auxSNP3Value = dev2Ctrls[iterId3 + i * numSNPs];

            aux = ctrls00[i] & auxSNP3Value;
            tripleCtrls[18] += _devpopcount(aux);

            aux = ctrls01[i] & auxSNP3Value;
            tripleCtrls[19] += _devpopcount(aux);

            aux = ctrls02[i] & auxSNP3Value;
            tripleCtrls[20] += _devpopcount(aux);

            aux = ctrls10[i] & auxSNP3Value;
            tripleCtrls[21] += _devpopcount(aux);

            aux = ctrls11[i] & auxSNP3Value;
            tripleCtrls[22] += _devpopcount(aux);

            aux = ctrls12[i] & auxSNP3Value;
            tripleCtrls[23] += _devpopcount(aux);

            aux = ctrls20[i] & auxSNP3Value;
            tripleCtrls[24] += _devpopcount(aux);

            aux = ctrls21[i] & auxSNP3Value;
            tripleCtrls[25] += _devpopcount(aux);

            aux = ctrls22[i] & auxSNP3Value;
            tripleCtrls[26] += _devpopcount(aux);
        }

        // Calculate the MI with the values of the contingency table
        float entX = 0.0;
        float entAll = 0.0;
        float pCase, pCtrl;

        for (int i = 0; i < 27; i++) {
            pCase = tripleCases[i] * _invInds;
            if (pCase != 0.0) {
                entAll -= pCase * log2(pCase);
            }

            pCtrl = tripleCtrls[i] * _invInds;
            if (pCtrl != 0.0) {
                entAll -= pCtrl * log2(pCtrl);
            }

            pCase += pCtrl;
            if (pCase != 0.0) {
                entX -= pCase * log2(pCase);
            }
        }

        // The result of the MI is now in entX
        entX += _entY - entAll;

#ifdef DEBUG
        printf("Thread %d in block %d: MI for triple (%u, %u, %u) is %f\n", threadIdx.x,
                blockIdx.x, myId1, myId2, iterId3, entX);
#endif

        // Now include the value in the output list if it is high enough
        // There are empty values in the array
        if (numEntriesWithMI < numOutputs) {
            myOutIds[numEntriesWithMI] = iterId3;
            myOutValues[numEntriesWithMI] = entX;

            // If this is the minimum value of the array
            if (entX < minMI) {
                minMI = entX;
                minMIPos = numEntriesWithMI;
            }

            numEntriesWithMI++;
        } else if (entX > minMI) { // The value must be inserted
            myOutIds[minMIPos] = iterId3;
            myOutValues[minMIPos] = entX;

            // Find the new minimum
            minMIPos = 0;
            minMI = myOutIds[0];
            for (int i = 1; i < numOutputs; i++) {
                if (myOutValues[i] < minMI) {
                    minMI = myOutValues[i];
                    minMIPos = i;
                }
            }
        }
    }

    // The thread has a list of numOutputs with the highest values
    // Complete the list just in case there are no so many values in total
    for (int i = numEntriesWithMI; i < numOutputs; i++) {
        myOutIds[i] = 0;
        myOutValues[i] = 0.0;
        minMI = 0.0;
        minMIPos = i;
    }

#ifdef DEBUG
    printf("Before reducing thread %d of block %d: %u (%f), %u (%f), %u (%f), %u (%f)\n",
                threadIdx.x, blockIdx.x,
                myOutIds[0], myOutValues[0], myOutIds[1], myOutValues[1],
                myOutIds[2], myOutValues[2], myOutIds[3], myOutValues[3]);
#endif

    float *remoteOutValues;
    uint32_t * remoteOutIds;

    // Perform the reduction of the lists of the block of threads
    // Each reduction obtains the numOutputs highest elements of two threads
    for (int stride = blockDim.x / 2; stride > 0; stride /= 2) {
        __syncthreads();
        if (threadIdx.x < stride) {
            // Each thread has its own minimum so we only need to compare the numOutputs values of the other thread
            remoteOutValues = &shMIValues[(threadIdx.x + stride) * numOutputs];
            remoteOutIds = &shMIId[(threadIdx.x + stride) * numOutputs];

            for (int i = 0; i < numOutputs; i++) {
                if (remoteOutValues[i] > minMI) { // The value must be inserted
                    myOutIds[minMIPos] = remoteOutIds[i];
                    myOutValues[minMIPos] = remoteOutValues[i];

                    // Find the new minimum
                    minMIPos = 0;
                    minMI = myOutIds[0];
                    for (int j = 1; j < numOutputs; j++) {
                        if (myOutValues[j] < minMI) {
                            minMI = myOutValues[j];
                            minMIPos = j;
                        }
                    }
                }
            }
#ifdef DEBUG
            printf("With stride %d thread %d of block %d: %u (%f), %u (%f), %u (%f), %u (%f)\n", stride,
                threadIdx.x, blockIdx.x,
                myOutIds[0], myOutValues[0], myOutIds[1], myOutValues[1],
                myOutIds[2], myOutValues[2], myOutIds[3], myOutValues[3]);
#endif
        }
    }

    // Save the output list for the block
    if (!threadIdx.x) {
        uint3 *blockOutIds = &devMiIds[blockIdx.x * numOutputs];
        float *blockOutValues = &devMIValues[blockIdx.x * numOutputs];

        for (int i = 0; i < numOutputs; i++) {
            blockOutIds[i].x = myId1;
            blockOutIds[i].y = myId2;
            blockOutIds[i].z = myOutIds[i];
            blockOutValues[i] = myOutValues[i];
        }
    }
}

static __global__ void _kernelTripleIG(uint64_t numPairs, uint32_t numSNPs,
                                       uint16_t numEntriesCases, uint16_t numEntriesCtrls, uint2 *devIds,
                                       uint32_t *dev0Cases, uint32_t *dev1Cases, uint32_t *dev2Cases,
                                       uint32_t *dev0Ctrls, uint32_t *dev1Ctrls, uint32_t *dev2Ctrls,
                                       GPUDoubleContTable *devDoubleTables, uint16_t numOutputs,
                                       float *devMIValues, uint3 *devMiIds) {

    extern __shared__ uint32_t
    sharedMem[];
    uint32_t * cases00 = sharedMem;
    uint32_t * cases01 = &sharedMem[numEntriesCases];
    uint32_t * cases02 = &sharedMem[2 * numEntriesCases];
    uint32_t * cases10 = &sharedMem[3 * numEntriesCases];
    uint32_t * cases11 = &sharedMem[4 * numEntriesCases];
    uint32_t * cases12 = &sharedMem[5 * numEntriesCases];
    uint32_t * cases20 = &sharedMem[6 * numEntriesCases];
    uint32_t * cases21 = &sharedMem[7 * numEntriesCases];
    uint32_t * cases22 = &sharedMem[8 * numEntriesCases];
    uint32_t * ctrls00 = &sharedMem[9 * numEntriesCases];
    uint32_t * ctrls01 = &sharedMem[9 * numEntriesCases + numEntriesCtrls];
    uint32_t * ctrls02 = &sharedMem[9 * numEntriesCases + 2 * numEntriesCtrls];
    uint32_t * ctrls10 = &sharedMem[9 * numEntriesCases + 3 * numEntriesCtrls];
    uint32_t * ctrls11 = &sharedMem[9 * numEntriesCases + 4 * numEntriesCtrls];
    uint32_t * ctrls12 = &sharedMem[9 * numEntriesCases + 5 * numEntriesCtrls];
    uint32_t * ctrls20 = &sharedMem[9 * numEntriesCases + 6 * numEntriesCtrls];
    uint32_t * ctrls21 = &sharedMem[9 * numEntriesCases + 7 * numEntriesCtrls];
    uint32_t * ctrls22 = &sharedMem[9 * numEntriesCases + 8 * numEntriesCtrls];
    uint32_t * shMIId = &sharedMem[9 * (numEntriesCases + numEntriesCtrls)];
    float *shMIValues = (float *) &sharedMem[9 * (numEntriesCases + numEntriesCtrls) + blockDim.x * numOutputs];

    uint32_t * myOutIds = &shMIId[threadIdx.x * numOutputs];
    float *myOutValues = &shMIValues[threadIdx.x * numOutputs];

    uint16_t numEntriesWithMI = 0;
    float minMI = _MAX_FLOAT;
    uint16_t minMIPos = 0;

    // Copy the information of the double contingency table to shared memory
    GPUDoubleContTable *table = &devDoubleTables[blockIdx.x];

    for (int iter = threadIdx.x; iter < numEntriesCases; iter += blockDim.x) {
        cases00[iter] = table->_cases00[iter];
        cases01[iter] = table->_cases01[iter];
        cases02[iter] = table->_cases02[iter];
        cases10[iter] = table->_cases10[iter];
        cases11[iter] = table->_cases11[iter];
        cases12[iter] = table->_cases12[iter];
        cases20[iter] = table->_cases20[iter];
        cases21[iter] = table->_cases21[iter];
        cases22[iter] = table->_cases22[iter];
    }
    for (int iter = threadIdx.x; iter < numEntriesCtrls; iter += blockDim.x) {
        ctrls00[iter] = table->_ctrls00[iter];
        ctrls01[iter] = table->_ctrls01[iter];
        ctrls02[iter] = table->_ctrls02[iter];
        ctrls10[iter] = table->_ctrls10[iter];
        ctrls11[iter] = table->_ctrls11[iter];
        ctrls12[iter] = table->_ctrls12[iter];
        ctrls20[iter] = table->_ctrls20[iter];
        ctrls21[iter] = table->_ctrls21[iter];
        ctrls22[iter] = table->_ctrls22[iter];
    }

    __syncthreads();

    uint16_t tripleCases[27];
    uint16_t tripleCtrls[27];

    uint32_t
    myId1 = devIds[blockIdx.x].x;
    uint32_t
    myId2 = devIds[blockIdx.x].y;
    uint32_t
    iterId3 = myId2 + threadIdx.x + 1;

    uint32_t
    aux;
    uint32_t
    auxSNP3Value;

    // Each thread computes several triples using the same pair
    for (; iterId3 < numSNPs; iterId3 += blockDim.x) {
        // Starts creating the contingency table
        for (int i = 0; i < 27; i++) {
            tripleCases[i] = 0;
            tripleCtrls[i] = 0;
        }

        for (int i = 0; i < numEntriesCases; i++) {
            auxSNP3Value = dev0Cases[iterId3 + i * numSNPs];

            aux = cases00[i] & auxSNP3Value;
            tripleCases[0] += _devpopcount(aux);

            aux = cases01[i] & auxSNP3Value;
            tripleCases[1] += _devpopcount(aux);

            aux = cases02[i] & auxSNP3Value;
            tripleCases[2] += _devpopcount(aux);

            aux = cases10[i] & auxSNP3Value;
            tripleCases[3] += _devpopcount(aux);

            aux = cases11[i] & auxSNP3Value;
            tripleCases[4] += _devpopcount(aux);

            aux = cases12[i] & auxSNP3Value;
            tripleCases[5] += _devpopcount(aux);

            aux = cases20[i] & auxSNP3Value;
            tripleCases[6] += _devpopcount(aux);

            aux = cases21[i] & auxSNP3Value;
            tripleCases[7] += _devpopcount(aux);

            aux = cases22[i] & auxSNP3Value;
            tripleCases[8] += _devpopcount(aux);

            auxSNP3Value = dev1Cases[iterId3 + i * numSNPs];

            aux = cases00[i] & auxSNP3Value;
            tripleCases[9] += _devpopcount(aux);

            aux = cases01[i] & auxSNP3Value;
            tripleCases[10] += _devpopcount(aux);

            aux = cases02[i] & auxSNP3Value;
            tripleCases[11] += _devpopcount(aux);

            aux = cases10[i] & auxSNP3Value;
            tripleCases[12] += _devpopcount(aux);

            aux = cases11[i] & auxSNP3Value;
            tripleCases[13] += _devpopcount(aux);

            aux = cases12[i] & auxSNP3Value;
            tripleCases[14] += _devpopcount(aux);

            aux = cases20[i] & auxSNP3Value;
            tripleCases[15] += _devpopcount(aux);

            aux = cases21[i] & auxSNP3Value;
            tripleCases[16] += _devpopcount(aux);

            aux = cases22[i] & auxSNP3Value;
            tripleCases[17] += _devpopcount(aux);

            auxSNP3Value = dev2Cases[iterId3 + i * numSNPs];

            aux = cases00[i] & auxSNP3Value;
            tripleCases[18] += _devpopcount(aux);

            aux = cases01[i] & auxSNP3Value;
            tripleCases[19] += _devpopcount(aux);

            aux = cases02[i] & auxSNP3Value;
            tripleCases[20] += _devpopcount(aux);

            aux = cases10[i] & auxSNP3Value;
            tripleCases[21] += _devpopcount(aux);

            aux = cases11[i] & auxSNP3Value;
            tripleCases[22] += _devpopcount(aux);

            aux = cases12[i] & auxSNP3Value;
            tripleCases[23] += _devpopcount(aux);

            aux = cases20[i] & auxSNP3Value;
            tripleCases[24] += _devpopcount(aux);

            aux = cases21[i] & auxSNP3Value;
            tripleCases[25] += _devpopcount(aux);

            aux = cases22[i] & auxSNP3Value;
            tripleCases[26] += _devpopcount(aux);
        }

        for (int i = 0; i < numEntriesCtrls; i++) {
            auxSNP3Value = dev0Ctrls[iterId3 + i * numSNPs];

            aux = ctrls00[i] & auxSNP3Value;
            tripleCtrls[0] += _devpopcount(aux);

            aux = ctrls01[i] & auxSNP3Value;
            tripleCtrls[1] += _devpopcount(aux);

            aux = ctrls02[i] & auxSNP3Value;
            tripleCtrls[2] += _devpopcount(aux);

            aux = ctrls10[i] & auxSNP3Value;
            tripleCtrls[3] += _devpopcount(aux);

            aux = ctrls11[i] & auxSNP3Value;
            tripleCtrls[4] += _devpopcount(aux);

            aux = ctrls12[i] & auxSNP3Value;
            tripleCtrls[5] += _devpopcount(aux);

            aux = ctrls20[i] & auxSNP3Value;
            tripleCtrls[6] += _devpopcount(aux);

            aux = ctrls21[i] & auxSNP3Value;
            tripleCtrls[7] += _devpopcount(aux);

            aux = ctrls22[i] & auxSNP3Value;
            tripleCtrls[8] += _devpopcount(aux);

            auxSNP3Value = dev1Ctrls[iterId3 + i * numSNPs];

            aux = ctrls00[i] & auxSNP3Value;
            tripleCtrls[9] += _devpopcount(aux);

            aux = ctrls01[i] & auxSNP3Value;
            tripleCtrls[10] += _devpopcount(aux);

            aux = ctrls02[i] & auxSNP3Value;
            tripleCtrls[11] += _devpopcount(aux);

            aux = ctrls10[i] & auxSNP3Value;
            tripleCtrls[12] += _devpopcount(aux);

            aux = ctrls11[i] & auxSNP3Value;
            tripleCtrls[13] += _devpopcount(aux);

            aux = ctrls12[i] & auxSNP3Value;
            tripleCtrls[14] += _devpopcount(aux);

            aux = ctrls20[i] & auxSNP3Value;
            tripleCtrls[15] += _devpopcount(aux);

            aux = ctrls21[i] & auxSNP3Value;
            tripleCtrls[16] += _devpopcount(aux);

            aux = ctrls22[i] & auxSNP3Value;
            tripleCtrls[17] += _devpopcount(aux);

            auxSNP3Value = dev2Ctrls[iterId3 + i * numSNPs];

            aux = ctrls00[i] & auxSNP3Value;
            tripleCtrls[18] += _devpopcount(aux);

            aux = ctrls01[i] & auxSNP3Value;
            tripleCtrls[19] += _devpopcount(aux);

            aux = ctrls02[i] & auxSNP3Value;
            tripleCtrls[20] += _devpopcount(aux);

            aux = ctrls10[i] & auxSNP3Value;
            tripleCtrls[21] += _devpopcount(aux);

            aux = ctrls11[i] & auxSNP3Value;
            tripleCtrls[22] += _devpopcount(aux);

            aux = ctrls12[i] & auxSNP3Value;
            tripleCtrls[23] += _devpopcount(aux);

            aux = ctrls20[i] & auxSNP3Value;
            tripleCtrls[24] += _devpopcount(aux);

            aux = ctrls21[i] & auxSNP3Value;
            tripleCtrls[25] += _devpopcount(aux);

            aux = ctrls22[i] & auxSNP3Value;
            tripleCtrls[26] += _devpopcount(aux);
        }

        float casex0 = 0.0;
        float casex1 = 0.0;
        float casex2 = 0.0;
        float casey0 = 0.0;
        float casey1 = 0.0;
        float casey2 = 0.0;
        float casez0 = 0.0;
        float casez1 = 0.0;
        float casez2 = 0.0;
        float casexy00 = 0.0;
        float casexy01 = 0.0;
        float casexy02 = 0.0;
        float casexy10 = 0.0;
        float casexy11 = 0.0;
        float casexy12 = 0.0;
        float casexy20 = 0.0;
        float casexy21 = 0.0;
        float casexy22 = 0.0;
        float casexz00 = 0.0;
        float casexz01 = 0.0;
        float casexz02 = 0.0;
        float casexz10 = 0.0;
        float casexz11 = 0.0;
        float casexz12 = 0.0;
        float casexz20 = 0.0;
        float casexz21 = 0.0;
        float casexz22 = 0.0;
        float caseyz00 = 0.0;
        float caseyz01 = 0.0;
        float caseyz02 = 0.0;
        float caseyz10 = 0.0;
        float caseyz11 = 0.0;
        float caseyz12 = 0.0;
        float caseyz20 = 0.0;
        float caseyz21 = 0.0;
        float caseyz22 = 0.0;
        float ctrlx0 = 0.0;
        float ctrlx1 = 0.0;
        float ctrlx2 = 0.0;
        float ctrly0 = 0.0;
        float ctrly1 = 0.0;
        float ctrly2 = 0.0;
        float ctrlz0 = 0.0;
        float ctrlz1 = 0.0;
        float ctrlz2 = 0.0;
        float ctrlxy00 = 0.0;
        float ctrlxy01 = 0.0;
        float ctrlxy02 = 0.0;
        float ctrlxy10 = 0.0;
        float ctrlxy11 = 0.0;
        float ctrlxy12 = 0.0;
        float ctrlxy20 = 0.0;
        float ctrlxy21 = 0.0;
        float ctrlxy22 = 0.0;
        float ctrlxz00 = 0.0;
        float ctrlxz01 = 0.0;
        float ctrlxz02 = 0.0;
        float ctrlxz10 = 0.0;
        float ctrlxz11 = 0.0;
        float ctrlxz12 = 0.0;
        float ctrlxz20 = 0.0;
        float ctrlxz21 = 0.0;
        float ctrlxz22 = 0.0;
        float ctrlyz00 = 0.0;
        float ctrlyz01 = 0.0;
        float ctrlyz02 = 0.0;
        float ctrlyz10 = 0.0;
        float ctrlyz11 = 0.0;
        float ctrlyz12 = 0.0;
        float ctrlyz20 = 0.0;
        float ctrlyz21 = 0.0;
        float ctrlyz22 = 0.0;
        float entAll = 0.0, entX = 0.0;
        float pCase, pCtrl, miXYZ;

        // Do it one by one in order to not repeat calculations or do extra if-else
        // Calculate the entropies
        // 000
        pCase = tripleCases[0] * _invInds;
        if (pCase != 0.0) {
            entAll -= pCase * log2(pCase);
        }
        casex0 += pCase;
        casey0 += pCase;
        casez0 += pCase;
        casexy00 += pCase;
        casexz00 += pCase;
        caseyz00 += pCase;

        pCtrl = tripleCtrls[0] * _invInds;
        if (pCtrl != 0.0) {
            entAll -= pCtrl * log2(pCtrl);
        }
        ctrlx0 += pCtrl;
        ctrly0 += pCtrl;
        ctrlz0 += pCtrl;
        ctrlxy00 += pCtrl;
        ctrlxz00 += pCtrl;
        ctrlyz00 += pCtrl;

        pCase += pCtrl;
        if (pCase != 0.0) {
            entX -= pCase * log2(pCase);
        }

        // 001
        pCase = tripleCases[1] * _invInds;
        if (pCase != 0.0) {
            entAll -= pCase * log2(pCase);
        }
        casex0 += pCase;
        casey0 += pCase;
        casez1 += pCase;
        casexy00 += pCase;
        casexz01 += pCase;
        caseyz01 += pCase;

        pCtrl = tripleCtrls[1] * _invInds;
        if (pCtrl != 0.0) {
            entAll -= pCtrl * log2(pCtrl);
        }
        ctrlx0 += pCtrl;
        ctrly0 += pCtrl;
        ctrlz1 += pCtrl;
        ctrlxy00 += pCtrl;
        ctrlxz01 += pCtrl;
        ctrlyz01 += pCtrl;

        pCase += pCtrl;
        if (pCase != 0.0) {
            entX -= pCase * log2(pCase);
        }

        // 002
        pCase = tripleCases[2] * _invInds;
        if (pCase != 0.0) {
            entAll -= pCase * log2(pCase);
        }
        casex0 += pCase;
        casey0 += pCase;
        casez2 += pCase;
        casexy00 += pCase;
        casexz02 += pCase;
        caseyz02 += pCase;

        pCtrl = tripleCtrls[2] * _invInds;
        if (pCtrl != 0.0) {
            entAll -= pCtrl * log2(pCtrl);
        }
        ctrlx0 += pCtrl;
        ctrly0 += pCtrl;
        ctrlz2 += pCtrl;
        ctrlxy00 += pCtrl;
        ctrlxz02 += pCtrl;
        ctrlyz02 += pCtrl;

        pCase += pCtrl;
        if (pCase != 0.0) {
            entX -= pCase * log2(pCase);
        }

        // 010
        pCase = tripleCases[3] * _invInds;
        if (pCase != 0.0) {
            entAll -= pCase * log2(pCase);
        }
        casex0 += pCase;
        casey1 += pCase;
        casez0 += pCase;
        casexy01 += pCase;
        casexz00 += pCase;
        caseyz10 += pCase;

        pCtrl = tripleCtrls[3] * _invInds;
        if (pCtrl != 0.0) {
            entAll -= pCtrl * log2(pCtrl);
        }
        ctrlx0 += pCtrl;
        ctrly1 += pCtrl;
        ctrlz0 += pCtrl;
        ctrlxy01 += pCtrl;
        ctrlxz00 += pCtrl;
        ctrlyz10 += pCtrl;

        pCase += pCtrl;
        if (pCase != 0.0) {
            entX -= pCase * log2(pCase);
        }

        // 011
        pCase = tripleCases[4] * _invInds;
        if (pCase != 0.0) {
            entAll -= pCase * log2(pCase);
        }
        casex0 += pCase;
        casey1 += pCase;
        casez1 += pCase;
        casexy01 += pCase;
        casexz01 += pCase;
        caseyz11 += pCase;

        pCtrl = tripleCtrls[4] * _invInds;
        if (pCtrl != 0.0) {
            entAll -= pCtrl * log2(pCtrl);
        }
        ctrlx0 += pCtrl;
        ctrly1 += pCtrl;
        ctrlz1 += pCtrl;
        ctrlxy01 += pCtrl;
        ctrlxz01 += pCtrl;
        ctrlyz11 += pCtrl;

        pCase += pCtrl;
        if (pCase != 0.0) {
            entX -= pCase * log2(pCase);
        }

        // 012
        pCase = tripleCases[5] * _invInds;
        if (pCase != 0.0) {
            entAll -= pCase * log2(pCase);
        }
        casex0 += pCase;
        casey1 += pCase;
        casez2 += pCase;
        casexy01 += pCase;
        casexz02 += pCase;
        caseyz12 += pCase;

        pCtrl = tripleCtrls[5] * _invInds;
        if (pCtrl != 0.0) {
            entAll -= pCtrl * log2(pCtrl);
        }
        ctrlx0 += pCtrl;
        ctrly1 += pCtrl;
        ctrlz2 += pCtrl;
        ctrlxy01 += pCtrl;
        ctrlxz02 += pCtrl;
        ctrlyz12 += pCtrl;

        pCase += pCtrl;
        if (pCase != 0.0) {
            entX -= pCase * log2(pCase);
        }

        // 020
        pCase = tripleCases[6] * _invInds;
        if (pCase != 0.0) {
            entAll -= pCase * log2(pCase);
        }
        casex0 += pCase;
        casey2 += pCase;
        casez0 += pCase;
        casexy02 += pCase;
        casexz00 += pCase;
        caseyz20 += pCase;

        pCtrl = tripleCtrls[6] * _invInds;
        if (pCtrl != 0.0) {
            entAll -= pCtrl * log2(pCtrl);
        }
        ctrlx0 += pCtrl;
        ctrly2 += pCtrl;
        ctrlz0 += pCtrl;
        ctrlxy02 += pCtrl;
        ctrlxz00 += pCtrl;
        ctrlyz20 += pCtrl;

        pCase += pCtrl;
        if (pCase != 0.0) {
            entX -= pCase * log2(pCase);
        }

        // 021
        pCase = tripleCases[7] * _invInds;
        if (pCase != 0.0) {
            entAll -= pCase * log2(pCase);
        }
        casex0 += pCase;
        casey2 += pCase;
        casez1 += pCase;
        casexy02 += pCase;
        casexz01 += pCase;
        caseyz21 += pCase;

        pCtrl = tripleCtrls[7] * _invInds;
        if (pCtrl != 0.0) {
            entAll -= pCtrl * log2(pCtrl);
        }
        ctrlx0 += pCtrl;
        ctrly2 += pCtrl;
        ctrlz1 += pCtrl;
        ctrlxy02 += pCtrl;
        ctrlxz01 += pCtrl;
        ctrlyz21 += pCtrl;

        pCase += pCtrl;
        if (pCase != 0.0) {
            entX -= pCase * log2(pCase);
        }

        // 022
        pCase = tripleCases[8] * _invInds;
        if (pCase != 0.0) {
            entAll -= pCase * log2(pCase);
        }
        casex0 += pCase;
        casey2 += pCase;
        casez2 += pCase;
        casexy02 += pCase;
        casexz02 += pCase;
        caseyz22 += pCase;

        pCtrl = tripleCtrls[8] * _invInds;
        if (pCtrl != 0.0) {
            entAll -= pCtrl * log2(pCtrl);
        }
        ctrlx0 += pCtrl;
        ctrly2 += pCtrl;
        ctrlz2 += pCtrl;
        ctrlxy02 += pCtrl;
        ctrlxz02 += pCtrl;
        ctrlyz22 += pCtrl;

        pCase += pCtrl;
        if (pCase != 0.0) {
            entX -= pCase * log2(pCase);
        }

        // 100
        pCase = tripleCases[9] * _invInds;
        if (pCase != 0.0) {
            entAll -= pCase * log2(pCase);
        }
        casex1 += pCase;
        casey0 += pCase;
        casez0 += pCase;
        casexy10 += pCase;
        casexz10 += pCase;
        caseyz00 += pCase;

        pCtrl = tripleCtrls[9] * _invInds;
        if (pCtrl != 0.0) {
            entAll -= pCtrl * log2(pCtrl);
        }
        ctrlx1 += pCtrl;
        ctrly0 += pCtrl;
        ctrlz0 += pCtrl;
        ctrlxy10 += pCtrl;
        ctrlxz10 += pCtrl;
        ctrlyz00 += pCtrl;

        pCase += pCtrl;
        if (pCase != 0.0) {
            entX -= pCase * log2(pCase);
        }

        // 101
        pCase = tripleCases[10] * _invInds;
        if (pCase != 0.0) {
            entAll -= pCase * log2(pCase);
        }
        casex1 += pCase;
        casey0 += pCase;
        casez1 += pCase;
        casexy10 += pCase;
        casexz11 += pCase;
        caseyz01 += pCase;

        pCtrl = tripleCtrls[10] * _invInds;
        if (pCtrl != 0.0) {
            entAll -= pCtrl * log2(pCtrl);
        }
        ctrlx1 += pCtrl;
        ctrly0 += pCtrl;
        ctrlz1 += pCtrl;
        ctrlxy10 += pCtrl;
        ctrlxz11 += pCtrl;
        ctrlyz01 += pCtrl;

        pCase += pCtrl;
        if (pCase != 0.0) {
            entX -= pCase * log2(pCase);
        }

        // 102
        pCase = tripleCases[11] * _invInds;
        if (pCase != 0.0) {
            entAll -= pCase * log2(pCase);
        }
        casex1 += pCase;
        casey0 += pCase;
        casez2 += pCase;
        casexy10 += pCase;
        casexz12 += pCase;
        caseyz02 += pCase;

        pCtrl = tripleCtrls[11] * _invInds;
        if (pCtrl != 0.0) {
            entAll -= pCtrl * log2(pCtrl);
        }
        ctrlx1 += pCtrl;
        ctrly0 += pCtrl;
        ctrlz2 += pCtrl;
        ctrlxy10 += pCtrl;
        ctrlxz12 += pCtrl;
        ctrlyz02 += pCtrl;

        pCase += pCtrl;
        if (pCase != 0.0) {
            entX -= pCase * log2(pCase);
        }

        // 110
        pCase = tripleCases[12] * _invInds;
        if (pCase != 0.0) {
            entAll -= pCase * log2(pCase);
        }
        casex1 += pCase;
        casey1 += pCase;
        casez0 += pCase;
        casexy11 += pCase;
        casexz10 += pCase;
        caseyz10 += pCase;

        pCtrl = tripleCtrls[12] * _invInds;
        if (pCtrl != 0.0) {
            entAll -= pCtrl * log2(pCtrl);
        }
        ctrlx1 += pCtrl;
        ctrly1 += pCtrl;
        ctrlz0 += pCtrl;
        ctrlxy11 += pCtrl;
        ctrlxz10 += pCtrl;
        ctrlyz10 += pCtrl;

        pCase += pCtrl;
        if (pCase != 0.0) {
            entX -= pCase * log2(pCase);
        }

        // 111
        pCase = tripleCases[13] * _invInds;
        if (pCase != 0.0) {
            entAll -= pCase * log2(pCase);
        }
        casex1 += pCase;
        casey1 += pCase;
        casez1 += pCase;
        casexy11 += pCase;
        casexz11 += pCase;
        caseyz11 += pCase;

        pCtrl = tripleCtrls[13] * _invInds;
        if (pCtrl != 0.0) {
            entAll -= pCtrl * log2(pCtrl);
        }
        ctrlx1 += pCtrl;
        ctrly1 += pCtrl;
        ctrlz1 += pCtrl;
        ctrlxy11 += pCtrl;
        ctrlxz11 += pCtrl;
        ctrlyz11 += pCtrl;

        pCase += pCtrl;
        if (pCase != 0.0) {
            entX -= pCase * log2(pCase);
        }

        // 112
        pCase = tripleCases[14] * _invInds;
        if (pCase != 0.0) {
            entAll -= pCase * log2(pCase);
        }
        casex1 += pCase;
        casey1 += pCase;
        casez2 += pCase;
        casexy11 += pCase;
        casexz12 += pCase;
        caseyz12 += pCase;

        pCtrl = tripleCtrls[14] * _invInds;
        if (pCtrl != 0.0) {
            entAll -= pCtrl * log2(pCtrl);
        }
        ctrlx1 += pCtrl;
        ctrly1 += pCtrl;
        ctrlz2 += pCtrl;
        ctrlxy11 += pCtrl;
        ctrlxz12 += pCtrl;
        ctrlyz12 += pCtrl;

        pCase += pCtrl;
        if (pCase != 0.0) {
            entX -= pCase * log2(pCase);
        }

        // 120
        pCase = tripleCases[15] * _invInds;
        if (pCase != 0.0) {
            entAll -= pCase * log2(pCase);
        }
        casex1 += pCase;
        casey2 += pCase;
        casez0 += pCase;
        casexy12 += pCase;
        casexz10 += pCase;
        caseyz20 += pCase;

        pCtrl = tripleCtrls[15] * _invInds;
        if (pCtrl != 0.0) {
            entAll -= pCtrl * log2(pCtrl);
        }
        ctrlx1 += pCtrl;
        ctrly2 += pCtrl;
        ctrlz0 += pCtrl;
        ctrlxy12 += pCtrl;
        ctrlxz10 += pCtrl;
        ctrlyz20 += pCtrl;

        pCase += pCtrl;
        if (pCase != 0.0) {
            entX -= pCase * log2(pCase);
        }

        // 121
        pCase = tripleCases[16] * _invInds;
        if (pCase != 0.0) {
            entAll -= pCase * log2(pCase);
        }
        casex1 += pCase;
        casey2 += pCase;
        casez1 += pCase;
        casexy12 += pCase;
        casexz11 += pCase;
        caseyz21 += pCase;

        pCtrl = tripleCtrls[16] * _invInds;
        if (pCtrl != 0.0) {
            entAll -= pCtrl * log2(pCtrl);
        }
        ctrlx1 += pCtrl;
        ctrly2 += pCtrl;
        ctrlz1 += pCtrl;
        ctrlxy12 += pCtrl;
        ctrlxz11 += pCtrl;
        ctrlyz21 += pCtrl;

        pCase += pCtrl;
        if (pCase != 0.0) {
            entX -= pCase * log2(pCase);
        }

        // 122
        pCase = tripleCases[17] * _invInds;
        if (pCase != 0.0) {
            entAll -= pCase * log2(pCase);
        }
        casex1 += pCase;
        casey2 += pCase;
        casez2 += pCase;
        casexy12 += pCase;
        casexz12 += pCase;
        caseyz22 += pCase;

        pCtrl = tripleCtrls[17] * _invInds;
        if (pCtrl != 0.0) {
            entAll -= pCtrl * log2(pCtrl);
        }
        ctrlx1 += pCtrl;
        ctrly2 += pCtrl;
        ctrlz2 += pCtrl;
        ctrlxy12 += pCtrl;
        ctrlxz12 += pCtrl;
        ctrlyz22 += pCtrl;

        pCase += pCtrl;
        if (pCase != 0.0) {
            entX -= pCase * log2(pCase);
        }

        // 200
        pCase = tripleCases[18] * _invInds;
        if (pCase != 0.0) {
            entAll -= pCase * log2(pCase);
        }
        casex2 += pCase;
        casey0 += pCase;
        casez0 += pCase;
        casexy20 += pCase;
        casexz20 += pCase;
        caseyz00 += pCase;

        pCtrl = tripleCtrls[18] * _invInds;
        if (pCtrl != 0.0) {
            entAll -= pCtrl * log2(pCtrl);
        }
        ctrlx2 += pCtrl;
        ctrly0 += pCtrl;
        ctrlz0 += pCtrl;
        ctrlxy20 += pCtrl;
        ctrlxz20 += pCtrl;
        ctrlyz00 += pCtrl;

        pCase += pCtrl;
        if (pCase != 0.0) {
            entX -= pCase * log2(pCase);
        }

        // 201
        pCase = tripleCases[19] * _invInds;
        if (pCase != 0.0) {
            entAll -= pCase * log2(pCase);
        }
        casex2 += pCase;
        casey0 += pCase;
        casez1 += pCase;
        casexy20 += pCase;
        casexz21 += pCase;
        caseyz01 += pCase;

        pCtrl = tripleCtrls[19] * _invInds;
        if (pCtrl != 0.0) {
            entAll -= pCtrl * log2(pCtrl);
        }
        ctrlx2 += pCtrl;
        ctrly0 += pCtrl;
        ctrlz1 += pCtrl;
        ctrlxy20 += pCtrl;
        ctrlxz21 += pCtrl;
        ctrlyz01 += pCtrl;

        pCase += pCtrl;
        if (pCase != 0.0) {
            entX -= pCase * log2(pCase);
        }

        // 202
        pCase = tripleCases[20] * _invInds;
        if (pCase != 0.0) {
            entAll -= pCase * log2(pCase);
        }
        casex2 += pCase;
        casey0 += pCase;
        casez2 += pCase;
        casexy20 += pCase;
        casexz22 += pCase;
        caseyz02 += pCase;

        pCtrl = tripleCtrls[20] * _invInds;
        if (pCtrl != 0.0) {
            entAll -= pCtrl * log2(pCtrl);
        }
        ctrlx2 += pCtrl;
        ctrly0 += pCtrl;
        ctrlz2 += pCtrl;
        ctrlxy20 += pCtrl;
        ctrlxz22 += pCtrl;
        ctrlyz02 += pCtrl;

        pCase += pCtrl;
        if (pCase != 0.0) {
            entX -= pCase * log2(pCase);
        }

        // 210
        pCase = tripleCases[21] * _invInds;
        if (pCase != 0.0) {
            entAll -= pCase * log2(pCase);
        }
        casex2 += pCase;
        casey1 += pCase;
        casez0 += pCase;
        casexy21 += pCase;
        casexz20 += pCase;
        caseyz10 += pCase;

        pCtrl = tripleCtrls[21] * _invInds;
        if (pCtrl != 0.0) {
            entAll -= pCtrl * log2(pCtrl);
        }
        ctrlx2 += pCtrl;
        ctrly1 += pCtrl;
        ctrlz0 += pCtrl;
        ctrlxy21 += pCtrl;
        ctrlxz20 += pCtrl;
        ctrlyz10 += pCtrl;

        pCase += pCtrl;
        if (pCase != 0.0) {
            entX -= pCase * log2(pCase);
        }

        // 211
        pCase = tripleCases[22] * _invInds;
        if (pCase != 0.0) {
            entAll -= pCase * log2(pCase);
        }
        casex2 += pCase;
        casey1 += pCase;
        casez1 += pCase;
        casexy21 += pCase;
        casexz21 += pCase;
        caseyz11 += pCase;

        pCtrl = tripleCtrls[22] * _invInds;
        if (pCtrl != 0.0) {
            entAll -= pCtrl * log2(pCtrl);
        }
        ctrlx2 += pCtrl;
        ctrly1 += pCtrl;
        ctrlz1 += pCtrl;
        ctrlxy21 += pCtrl;
        ctrlxz21 += pCtrl;
        ctrlyz11 += pCtrl;

        pCase += pCtrl;
        if (pCase != 0.0) {
            entX -= pCase * log2(pCase);
        }

        // 212
        pCase = tripleCases[23] * _invInds;
        if (pCase != 0.0) {
            entAll -= pCase * log2(pCase);
        }
        casex2 += pCase;
        casey1 += pCase;
        casez2 += pCase;
        casexy21 += pCase;
        casexz22 += pCase;
        caseyz12 += pCase;

        pCtrl = tripleCtrls[23] * _invInds;
        if (pCtrl != 0.0) {
            entAll -= pCtrl * log2(pCtrl);
        }
        ctrlx2 += pCtrl;
        ctrly1 += pCtrl;
        ctrlz2 += pCtrl;
        ctrlxy21 += pCtrl;
        ctrlxz22 += pCtrl;
        ctrlyz12 += pCtrl;

        pCase += pCtrl;
        if (pCase != 0.0) {
            entX -= pCase * log2(pCase);
        }

        // 220
        pCase = tripleCases[24] * _invInds;
        if (pCase != 0.0) {
            entAll -= pCase * log2(pCase);
        }
        casex2 += pCase;
        casey2 += pCase;
        casez0 += pCase;
        casexy22 += pCase;
        casexz20 += pCase;
        caseyz20 += pCase;

        pCtrl = tripleCtrls[24] * _invInds;
        if (pCtrl != 0.0) {
            entAll -= pCtrl * log2(pCtrl);
        }
        ctrlx2 += pCtrl;
        ctrly2 += pCtrl;
        ctrlz0 += pCtrl;
        ctrlxy22 += pCtrl;
        ctrlxz20 += pCtrl;
        ctrlyz20 += pCtrl;

        pCase += pCtrl;
        if (pCase != 0.0) {
            entX -= pCase * log2(pCase);
        }

        // 221
        pCase = tripleCases[25] * _invInds;
        if (pCase != 0.0) {
            entAll -= pCase * log2(pCase);
        }
        casex2 += pCase;
        casey2 += pCase;
        casez1 += pCase;
        casexy22 += pCase;
        casexz21 += pCase;
        caseyz21 += pCase;

        pCtrl = tripleCtrls[25] * _invInds;
        if (pCtrl != 0.0) {
            entAll -= pCtrl * log2(pCtrl);
        }
        ctrlx2 += pCtrl;
        ctrly2 += pCtrl;
        ctrlz1 += pCtrl;
        ctrlxy22 += pCtrl;
        ctrlxz21 += pCtrl;
        ctrlyz21 += pCtrl;

        pCase += pCtrl;
        if (pCase != 0.0) {
            entX -= pCase * log2(pCase);
        }

        // 222
        pCase = tripleCases[26] * _invInds;
        if (pCase != 0.0) {
            entAll -= pCase * log2(pCase);
        }
        casex2 += pCase;
        casey2 += pCase;
        casez2 += pCase;
        casexy22 += pCase;
        casexz22 += pCase;
        caseyz22 += pCase;

        pCtrl = tripleCtrls[26] * _invInds;
        if (pCtrl != 0.0) {
            entAll -= pCtrl * log2(pCtrl);
        }
        ctrlx2 += pCtrl;
        ctrly2 += pCtrl;
        ctrlz2 += pCtrl;
        ctrlxy22 += pCtrl;
        ctrlxz22 += pCtrl;
        ctrlyz22 += pCtrl;

        pCase += pCtrl;
        if (pCase != 0.0) {
            entX -= pCase * log2(pCase);
        }


        float miX = _entY;
        if (casex0 != 0.0) {
            miX += casex0 * log2(casex0);
        }
        if (casex1 != 0.0) {
            miX += casex1 * log2(casex1);
        }
        if (casex2 != 0.0) {
            miX += casex2 * log2(casex2);
        }
        if (ctrlx0 != 0.0) {
            miX += ctrlx0 * log2(ctrlx0);
        }
        if (ctrlx1 != 0.0) {
            miX += ctrlx1 * log2(ctrlx1);
        }
        if (ctrlx2 != 0.0) {
            miX += ctrlx2 * log2(ctrlx2);
        }


        casex0 += ctrlx0;
        if (casex0 != 0.0) {
            miX -= casex0 * log2(casex0);
        }
        casex1 += ctrlx1;
        if (casex1 != 0.0) {
            miX -= casex1 * log2(casex1);
        }
        casex2 += ctrlx2;
        if (casex2 != 0.0) {
            miX -= casex2 * log2(casex2);
        }


        float miY = _entY;
        if (casey0 != 0.0) {
            miY += casey0 * log2(casey0);
        }
        if (casey1 != 0.0) {
            miY += casey1 * log2(casey1);
        }
        if (casey2 != 0.0) {
            miY += casey2 * log2(casey2);
        }
        if (ctrly0 != 0.0) {
            miY += ctrly0 * log2(ctrly0);
        }
        if (ctrly1 != 0.0) {
            miY += ctrly1 * log2(ctrly1);
        }
        if (ctrly2 != 0.0) {
            miY += ctrly2 * log2(ctrly2);
        }


        casey0 += ctrly0;
        if (casey0 != 0.0) {
            miY -= casey0 * log2(casey0);
        }
        casey1 += ctrly1;
        if (casey1 != 0.0) {
            miY -= casey1 * log2(casey1);
        }
        casey2 += ctrly2;
        if (casey2 != 0.0) {
            miY -= casey2 * log2(casey2);
        }


        float miZ = _entY;
        if (casez0 != 0.0) {
            miZ += casez0 * log2(casez0);
        }
        if (casez1 != 0.0) {
            miZ += casez1 * log2(casez1);
        }
        if (casez2 != 0.0) {
            miZ += casez2 * log2(casez2);
        }
        if (ctrlz0 != 0.0) {
            miZ += ctrlz0 * log2(ctrlz0);
        }
        if (ctrlz1 != 0.0) {
            miZ += ctrlz1 * log2(ctrlz1);
        }
        if (ctrlz2 != 0.0) {
            miZ += ctrlz2 * log2(ctrlz2);
        }


        casez0 += ctrlz0;
        if (casez0 != 0.0) {
            miZ -= casez0 * log2(casez0);
        }
        casez1 += ctrlz1;
        if (casez1 != 0.0) {
            miZ -= casez1 * log2(casez1);
        }
        casez2 += ctrlz2;
        if (casez2 != 0.0) {
            miZ -= casez2 * log2(casez2);
        }


        float igXY = _entY;
        if (casexy00 != 0.0) {
            igXY += casexy00 * log2(casexy00);
        }
        if (casexy01 != 0.0) {
            igXY += casexy01 * log2(casexy01);
        }
        if (casexy02 != 0.0) {
            igXY += casexy02 * log2(casexy02);
        }
        if (casexy10 != 0.0) {
            igXY += casexy10 * log2(casexy10);
        }
        if (casexy11 != 0.0) {
            igXY += casexy11 * log2(casexy11);
        }
        if (casexy12 != 0.0) {
            igXY += casexy12 * log2(casexy12);
        }
        if (casexy20 != 0.0) {
            igXY += casexy20 * log2(casexy20);
        }
        if (casexy21 != 0.0) {
            igXY += casexy21 * log2(casexy21);
        }
        if (casexy22 != 0.0) {
            igXY += casexy22 * log2(casexy22);
        }
        if (ctrlxy00 != 0.0) {
            igXY += ctrlxy00 * log2(ctrlxy00);
        }
        if (ctrlxy01 != 0.0) {
            igXY += ctrlxy01 * log2(ctrlxy01);
        }
        if (ctrlxy02 != 0.0) {
            igXY += ctrlxy02 * log2(ctrlxy02);
        }
        if (ctrlxy10 != 0.0) {
            igXY += ctrlxy10 * log2(ctrlxy10);
        }
        if (ctrlxy11 != 0.0) {
            igXY += ctrlxy11 * log2(ctrlxy11);
        }
        if (ctrlxy12 != 0.0) {
            igXY += ctrlxy12 * log2(ctrlxy12);
        }
        if (ctrlxy20 != 0.0) {
            igXY += ctrlxy20 * log2(ctrlxy20);
        }
        if (ctrlxy21 != 0.0) {
            igXY += ctrlxy21 * log2(ctrlxy21);
        }
        if (ctrlxy22 != 0.0) {
            igXY += ctrlxy22 * log2(ctrlxy22);
        }


        casexy00 += ctrlxy00;
        if (casexy00 != 0.0) {
            igXY -= casexy00 * log2(casexy00);
        }
        casexy01 += ctrlxy01;
        if (casexy01 != 0.0) {
            igXY -= casexy01 * log2(casexy01);
        }
        casexy02 += ctrlxy02;
        if (casexy02 != 0.0) {
            igXY -= casexy02 * log2(casexy02);
        }
        casexy10 += ctrlxy10;
        if (casexy10 != 0.0) {
            igXY -= casexy10 * log2(casexy10);
        }
        casexy11 += ctrlxy11;
        if (casexy11 != 0.0) {
            igXY -= casexy11 * log2(casexy11);
        }
        casexy12 += ctrlxy12;
        if (casexy12 != 0.0) {
            igXY -= casexy12 * log2(casexy12);
        }
        casexy20 += ctrlxy20;
        if (casexy20 != 0.0) {
            igXY -= casexy20 * log2(casexy20);
        }
        casexy21 += ctrlxy21;
        if (casexy21 != 0.0) {
            igXY -= casexy21 * log2(casexy21);
        }
        casexy22 += ctrlxy22;
        if (casexy22 != 0.0) {
            igXY -= casexy22 * log2(casexy22);
        }


        igXY -= miX + miY;
        if (igXY < 0.0) {
            igXY = 0.0;
        }


        float igXZ = _entY;
        if (casexz00 != 0.0) {
            igXZ += casexz00 * log2(casexz00);
        }
        if (casexz01 != 0.0) {
            igXZ += casexz01 * log2(casexz01);
        }
        if (casexz02 != 0.0) {
            igXZ += casexz02 * log2(casexz02);
        }
        if (casexz10 != 0.0) {
            igXZ += casexz10 * log2(casexz10);
        }
        if (casexz11 != 0.0) {
            igXZ += casexz11 * log2(casexz11);
        }
        if (casexz12 != 0.0) {
            igXZ += casexz12 * log2(casexz12);
        }
        if (casexz20 != 0.0) {
            igXZ += casexz20 * log2(casexz20);
        }
        if (casexz21 != 0.0) {
            igXZ += casexz21 * log2(casexz21);
        }
        if (casexz22 != 0.0) {
            igXZ += casexz22 * log2(casexz22);
        }
        if (ctrlxz00 != 0.0) {
            igXZ += ctrlxz00 * log2(ctrlxz00);
        }
        if (ctrlxz01 != 0.0) {
            igXZ += ctrlxz01 * log2(ctrlxz01);
        }
        if (ctrlxz02 != 0.0) {
            igXZ += ctrlxz02 * log2(ctrlxz02);
        }
        if (ctrlxz10 != 0.0) {
            igXZ += ctrlxz10 * log2(ctrlxz10);
        }
        if (ctrlxz11 != 0.0) {
            igXZ += ctrlxz11 * log2(ctrlxz11);
        }
        if (ctrlxz12 != 0.0) {
            igXZ += ctrlxz12 * log2(ctrlxz12);
        }
        if (ctrlxz20 != 0.0) {
            igXZ += ctrlxz20 * log2(ctrlxz20);
        }
        if (ctrlxz21 != 0.0) {
            igXZ += ctrlxz21 * log2(ctrlxz21);
        }
        if (ctrlxz22 != 0.0) {
            igXZ += ctrlxz22 * log2(ctrlxy22);
        }

        casexz00 += ctrlxz00;
        if (casexz00 != 0.0) {
            igXZ -= casexz00 * log2(casexz00);
        }
        casexz01 += ctrlxz01;
        if (casexz01 != 0.0) {
            igXZ -= casexz01 * log2(casexz01);
        }
        casexz02 += ctrlxz02;
        if (casexz02 != 0.0) {
            igXZ -= casexz02 * log2(casexz02);
        }
        casexz10 += ctrlxz10;
        if (casexz10 != 0.0) {
            igXZ -= casexz10 * log2(casexz10);
        }
        casexz11 += ctrlxz11;
        if (casexz11 != 0.0) {
            igXZ -= casexz11 * log2(casexz11);
        }
        casexz12 += ctrlxz12;
        if (casexz12 != 0.0) {
            igXZ -= casexz12 * log2(casexz12);
        }
        casexz20 += ctrlxz20;
        if (casexz20 != 0.0) {
            igXZ -= casexz20 * log2(casexz20);
        }
        casexz21 += ctrlxz21;
        if (casexz21 != 0.0) {
            igXZ -= casexz21 * log2(casexz21);
        }
        casexz22 += ctrlxz22;
        if (casexz22 != 0.0) {
            igXZ -= casexz22 * log2(casexz22);
        }

        igXZ -= miX + miZ;
        if (igXZ < 0.0) {
            igXZ = 0.0;
        }

        float igYZ = _entY;
        if (caseyz00 != 0.0) {
            igYZ += caseyz00 * log2(caseyz00);
        }
        if (caseyz01 != 0.0) {
            igYZ += caseyz01 * log2(caseyz01);
        }
        if (caseyz02 != 0.0) {
            igYZ += caseyz02 * log2(caseyz02);
        }
        if (caseyz10 != 0.0) {
            igYZ += caseyz10 * log2(caseyz10);
        }
        if (caseyz11 != 0.0) {
            igYZ += caseyz11 * log2(caseyz11);
        }
        if (caseyz12 != 0.0) {
            igYZ += caseyz12 * log2(caseyz12);
        }
        if (caseyz20 != 0.0) {
            igYZ += caseyz20 * log2(caseyz20);
        }
        if (caseyz21 != 0.0) {
            igYZ += caseyz21 * log2(caseyz21);
        }
        if (caseyz22 != 0.0) {
            igYZ += caseyz22 * log2(caseyz22);
        }
        if (ctrlyz00 != 0.0) {
            igYZ += ctrlyz00 * log2(ctrlyz00);
        }
        if (ctrlyz01 != 0.0) {
            igYZ += ctrlyz01 * log2(ctrlyz01);
        }
        if (ctrlyz02 != 0.0) {
            igYZ += ctrlyz02 * log2(ctrlyz02);
        }
        if (ctrlyz10 != 0.0) {
            igYZ += ctrlyz10 * log2(ctrlyz10);
        }
        if (ctrlyz11 != 0.0) {
            igYZ += ctrlyz11 * log2(ctrlyz11);
        }
        if (ctrlyz12 != 0.0) {
            igYZ += ctrlyz12 * log2(ctrlyz12);
        }
        if (ctrlyz20 != 0.0) {
            igYZ += ctrlyz20 * log2(ctrlyz20);
        }
        if (ctrlyz21 != 0.0) {
            igYZ += ctrlyz21 * log2(ctrlyz21);
        }
        if (ctrlyz22 != 0.0) {
            igYZ += ctrlyz22 * log2(ctrlyz22);
        }

        caseyz00 += ctrlyz00;
        if (caseyz00 != 0.0) {
            igYZ -= caseyz00 * log2(caseyz00);
        }
        caseyz01 += ctrlyz01;
        if (caseyz01 != 0.0) {
            igYZ -= caseyz01 * log2(caseyz01);
        }
        caseyz02 += ctrlyz02;
        if (caseyz02 != 0.0) {
            igYZ -= caseyz02 * log2(caseyz02);
        }
        caseyz10 += ctrlyz10;
        if (caseyz10 != 0.0) {
            igYZ -= caseyz10 * log2(caseyz10);
        }
        caseyz11 += ctrlyz11;
        if (caseyz11 != 0.0) {
            igYZ -= caseyz11 * log2(caseyz11);
        }
        caseyz12 += ctrlyz12;
        if (caseyz12 != 0.0) {
            igYZ -= caseyz12 * log2(caseyz12);
        }
        caseyz20 += ctrlyz20;
        if (caseyz20 != 0.0) {
            igYZ -= caseyz20 * log2(caseyz20);
        }
        caseyz21 += ctrlyz21;
        if (caseyz21 != 0.0) {
            igYZ -= caseyz21 * log2(caseyz21);
        }
        caseyz22 += ctrlyz22;
        if (caseyz22 != 0.0) {
            igYZ -= caseyz22 * log2(caseyz22);
        }

        igYZ -= miZ + miY;
        if (igYZ < 0.0) {
            igYZ = 0.0;
        }

        miXYZ = _entY + entX - entAll;

        //printf("Thread %d in block %d: miXYZ for triple (%u, %u, %u) is %f (%f+%f-%f)\n", threadIdx.x,
        //	blockIdx.x, myId1, myId2, iterId3, miXYZ, _entY, entX, entAll);

        entX = miXYZ - igXY - igXZ - igYZ - miX - miY - miZ;

#ifdef DEBUG
        printf("Thread %d in block %d: IG for triple (%u, %u, %u) is %f (%f-%f-%f-%f-%f-%f-%f)\n", threadIdx.x,
                blockIdx.x, myId1, myId2, iterId3, entX, miXYZ, igXY, igXZ, miX, miY, miZ);
#endif


        // Now include the value in the output list if it is high enough
        // There are empty values in the array
        if (numEntriesWithMI < numOutputs) {
            myOutIds[numEntriesWithMI] = iterId3;
            myOutValues[numEntriesWithMI] = entX;

            // If this is the minimum value of the array
            if (entX < minMI) {
                minMI = entX;
                minMIPos = numEntriesWithMI;
            }

            numEntriesWithMI++;
        } else if (entX > minMI) { // The value must be inserted
            myOutIds[minMIPos] = iterId3;
            myOutValues[minMIPos] = entX;

            // Find the new minimum
            minMIPos = 0;
            minMI = myOutIds[0];
            for (int i = 1; i < numOutputs; i++) {
                if (myOutValues[i] < minMI) {
                    minMI = myOutValues[i];
                    minMIPos = i;
                }
            }
        }
    }

    // The thread has a list of numOutputs with the highest values
    // Complete the list just in case there are no so many values in total
    for (int i = numEntriesWithMI; i < numOutputs; i++) {
        myOutIds[i] = 0;
        myOutValues[i] = 0.0;
        minMI = 0.0;
        minMIPos = i;
    }

#ifdef DEBUG
    printf("Before reducing thread %d of block %d: %u (%f), %u (%f), %u (%f), %u (%f)\n",
                threadIdx.x, blockIdx.x,
                myOutIds[0], myOutValues[0], myOutIds[1], myOutValues[1],
                myOutIds[2], myOutValues[2], myOutIds[3], myOutValues[3]);
#endif

    float *remoteOutValues;
    uint32_t * remoteOutIds;

    // Perform the reduction of the lists of the block of threads
    // Each reduction obtains the numOutputs highest elements of two threads
    for (int stride = blockDim.x / 2; stride > 0; stride /= 2) {
        __syncthreads();
        if (threadIdx.x < stride) {
            // Each thread has its own minimum so we only need to compare the numOutputs values of the other thread
            remoteOutValues = &shMIValues[(threadIdx.x + stride) * numOutputs];
            remoteOutIds = &shMIId[(threadIdx.x + stride) * numOutputs];

            for (int i = 0; i < numOutputs; i++) {
                if (remoteOutValues[i] > minMI) { // The value must be inserted
                    myOutIds[minMIPos] = remoteOutIds[i];
                    myOutValues[minMIPos] = remoteOutValues[i];

                    // Find the new minimum
                    minMIPos = 0;
                    minMI = myOutIds[0];
                    for (int j = 1; j < numOutputs; j++) {
                        if (myOutValues[j] < minMI) {
                            minMI = myOutValues[j];
                            minMIPos = j;
                        }
                    }
                }
            }
#ifdef DEBUG
            printf("With stride %d thread %d of block %d: %u (%f), %u (%f), %u (%f), %u (%f)\n", stride,
                threadIdx.x, blockIdx.x,
                myOutIds[0], myOutValues[0], myOutIds[1], myOutValues[1],
                myOutIds[2], myOutValues[2], myOutIds[3], myOutValues[3]);
#endif
        }
    }

    // Save the output list for the block
    if (!threadIdx.x) {
        uint3 *blockOutIds = &devMiIds[blockIdx.x * numOutputs];
        float *blockOutValues = &devMIValues[blockIdx.x * numOutputs];

        for (int i = 0; i < numOutputs; i++) {
            blockOutIds[i].x = myId1;
            blockOutIds[i].y = myId2;
            blockOutIds[i].z = myOutIds[i];
            blockOutValues[i] = myOutValues[i];
        }
    }
}

MutualInformation::MutualInformation(bool isMI, uint32_t numSNPs, uint16_t numCases, uint16_t numCtrls,
                                     std::vector<std::vector<uint32_t> *> cases,
                                     std::vector<std::vector<uint32_t> *> ctrls) :
        _isMI(isMI),
        _numSNPs(numSNPs),
        _numEntriesCase(numCases / 32 + ((numCases % 32) > 0)),
        _numEntriesCtrl(numCtrls / 32 + ((numCtrls % 32) > 0)) {
    // Allocate the arrays
    if (cudaSuccess != cudaMalloc(&_dev0Cases, _numEntriesCase * _numSNPs * sizeof(uint32_t)))
        throw CUDAError();
    if (cudaSuccess != cudaMalloc(&_dev1Cases, _numEntriesCase * _numSNPs * sizeof(uint32_t)))
        throw CUDAError();
    if (cudaSuccess != cudaMalloc(&_dev2Cases, _numEntriesCase * _numSNPs * sizeof(uint32_t)))
        throw CUDAError();
    if (cudaSuccess != cudaMalloc(&_dev0Ctrls, _numEntriesCtrl * _numSNPs * sizeof(uint32_t)))
        throw CUDAError();
    if (cudaSuccess != cudaMalloc(&_dev1Ctrls, _numEntriesCtrl * _numSNPs * sizeof(uint32_t)))
        throw CUDAError();
    if (cudaSuccess != cudaMalloc(&_dev2Ctrls, _numEntriesCtrl * _numSNPs * sizeof(uint32_t)))
        throw CUDAError();

    // All the entries are cyclicly ordered by SNPs
    if (cudaSuccess !=
        cudaMemcpy(_dev0Cases, &cases[0][0][0], _numSNPs * _numEntriesCase * sizeof(uint32_t), cudaMemcpyHostToDevice))
        throw CUDAError();
    if (cudaSuccess !=
        cudaMemcpy(_dev1Cases, &cases[0][1][0], _numSNPs * _numEntriesCase * sizeof(uint32_t), cudaMemcpyHostToDevice))
        throw CUDAError();
    if (cudaSuccess !=
        cudaMemcpy(_dev2Cases, &cases[0][2][0], _numSNPs * _numEntriesCase * sizeof(uint32_t), cudaMemcpyHostToDevice))
        throw CUDAError();
    if (cudaSuccess !=
        cudaMemcpy(_dev0Ctrls, &ctrls[0][0][0], _numSNPs * _numEntriesCtrl * sizeof(uint32_t), cudaMemcpyHostToDevice))
        throw CUDAError();
    if (cudaSuccess !=
        cudaMemcpy(_dev1Ctrls, &ctrls[0][1][0], _numSNPs * _numEntriesCtrl * sizeof(uint32_t), cudaMemcpyHostToDevice))
        throw CUDAError();
    if (cudaSuccess !=
        cudaMemcpy(_dev2Ctrls, &ctrls[0][2][0], _numSNPs * _numEntriesCtrl * sizeof(uint32_t), cudaMemcpyHostToDevice))
        throw CUDAError();

    // Increase the size of the L1 compared to shared memory
    if (cudaSuccess != cudaFuncSetCacheConfig(_kernelDoubleTable, cudaFuncCachePreferL1))
        throw CUDAError();

    // Increase the size of the shared memory compared to L1
    if (cudaSuccess != cudaFuncSetCacheConfig(_kernelTripleMI, cudaFuncCachePreferShared))
        throw CUDAError();

    // Increase the size of the shared memory compared to L1
    if (cudaSuccess != cudaFuncSetCacheConfig(_kernelTripleIG, cudaFuncCachePreferShared))
        throw CUDAError();

    float invInds = 1.0 / (numCases + numCtrls);
    float p = numCases * invInds;
    float entY = (-1.0) * p * log2(p);

    p = numCtrls * invInds;
    entY -= p * log2(p);

    if (cudaSuccess != cudaMemcpyToSymbol(_invInds, &invInds, sizeof(float), 0, cudaMemcpyHostToDevice))
        throw CUDAError();

    if (cudaSuccess != cudaMemcpyToSymbol(_entY, &entY, sizeof(float), 0, cudaMemcpyHostToDevice))
        throw CUDAError();

    float maxFL = FLT_MAX;
    if (cudaSuccess != cudaMemcpyToSymbol(_MAX_FLOAT, &maxFL, sizeof(float), 0, cudaMemcpyHostToDevice))
        throw CUDAError();
}

MutualInformation::~MutualInformation() {
    if (cudaSuccess != cudaFree(_dev0Cases))
        throw CUDAError();
    if (cudaSuccess != cudaFree(_dev1Cases))
        throw CUDAError();
    if (cudaSuccess != cudaFree(_dev2Cases))
        throw CUDAError();
    if (cudaSuccess != cudaFree(_dev0Ctrls))
        throw CUDAError();
    if (cudaSuccess != cudaFree(_dev1Ctrls))
        throw CUDAError();
    if (cudaSuccess != cudaFree(_dev2Ctrls))
        throw CUDAError();
}

long MutualInformation::compute(const std::vector<uint2> &pairs, uint16_t num_outputs, Position *output) {
    constexpr size_t block_size = 5000;

    if (cudaSuccess != cudaMalloc(&_devIds, block_size * sizeof(uint2)))
        throw CUDAError();

    // Auxiliary array for the contingency tables between the two kernels
    GPUDoubleContTable *_tables = new GPUDoubleContTable[block_size];
    for (int i = 0; i < block_size; i++) {
        _tables[i].initialize(_numEntriesCase, _numEntriesCtrl);
    }

    GPUDoubleContTable *_devDoubleTables;
    if (cudaSuccess != cudaMalloc(&_devDoubleTables, block_size * sizeof(GPUDoubleContTable)))
        throw CUDAError();
    if (cudaSuccess !=
        cudaMemcpy(_devDoubleTables, _tables, block_size * sizeof(GPUDoubleContTable),
                   cudaMemcpyHostToDevice))
        throw CUDAError();

    // Auxiliary array to store the MI values of each block
    float *_devMIValues;
    if (cudaSuccess != cudaMalloc(&_devMIValues, block_size * num_outputs * sizeof(float)))
        throw CUDAError();
    float *_hostMIValues = new float[block_size * num_outputs];

    // Auxiliary arrays to store the ids that are in the list of MIs
    uint3 *_devMiIds;
    if (cudaSuccess != cudaMalloc(&_devMiIds, block_size * num_outputs * sizeof(uint3)))
        throw CUDAError();
    uint3 *_hostMiIds = new uint3[block_size * num_outputs];

    // The minimum value in the array
    float minMI = FLT_MAX;
    // The position of the minimum value
    uint16_t minMIPos = 0;
    // Number of entries of the array full
    uint16_t numEntriesWithMI = 0;

    for (unsigned long i = 0; i < pairs.size(); i += block_size) {
        const unsigned long num_pairs = pairs.size() - i < block_size ? pairs.size() - i : block_size;

        if (cudaSuccess != cudaMemcpy(_devIds, &pairs.at(0) + i, num_pairs * sizeof(uint2), cudaMemcpyHostToDevice))
            throw CUDAError();

        uint32_t
        nblocks = (num_pairs + NUM_TH_PER_BLOCK - 1) / NUM_TH_PER_BLOCK;
        dim3 gridDouble(nblocks, 1);
        dim3 blocksDouble(NUM_TH_PER_BLOCK, 1);

        // Starts computing the double contingency tables
        _kernelDoubleTable << < gridDouble, blocksDouble, 0 >> > (num_pairs, _numSNPs,
                _numEntriesCase, _numEntriesCtrl, _devIds,
                _dev0Cases, _dev1Cases, _dev2Cases, _dev0Ctrls, _dev1Ctrls, _dev2Ctrls,
                _devDoubleTables);

        // Now you need to calculate the MI for each triple
        // Each block performs all the triples for one pair
        // The necessary shared memory is to store the double contingency table of the matrix
        dim3 gridMI(num_pairs, 1);
        dim3 blockMI(NUM_TH_PER_BLOCK, 1);
        uint32_t
        sharedSize = 9 * (_numEntriesCase + _numEntriesCtrl) * sizeof(uint32_t);
        sharedSize += num_outputs * NUM_TH_PER_BLOCK * (sizeof(float) + sizeof(uint32_t));

        if (_isMI) {
            _kernelTripleMI << < gridMI, blockMI, sharedSize >> >
                                                  (num_pairs, _numSNPs, _numEntriesCase, _numEntriesCtrl,
                                                          _devIds, _dev0Cases, _dev1Cases, _dev2Cases, _dev0Ctrls, _dev1Ctrls, _dev2Ctrls,
                                                          _devDoubleTables, num_outputs, _devMIValues, _devMiIds);
        } else {
            _kernelTripleIG << < gridMI, blockMI, sharedSize >> >
                                                  (num_pairs, _numSNPs, _numEntriesCase, _numEntriesCtrl,
                                                          _devIds, _dev0Cases, _dev1Cases, _dev2Cases, _dev0Ctrls, _dev1Ctrls, _dev2Ctrls,
                                                          _devDoubleTables, num_outputs, _devMIValues, _devMiIds);
        }

        if (cudaSuccess != cudaMemcpy(_hostMIValues, _devMIValues, num_pairs * num_outputs * sizeof(float),
                                      cudaMemcpyDeviceToHost))
            throw CUDAError();
        if (cudaSuccess != cudaMemcpy(_hostMiIds, _devMiIds, num_pairs * num_outputs * sizeof(uint3),
                                      cudaMemcpyDeviceToHost))
            throw CUDAError();

        _findNHighestMI(_hostMiIds, _hostMIValues, num_pairs * num_outputs, minMI, minMIPos, numEntriesWithMI,
                        num_outputs, output);
    }

    if (cudaSuccess != cudaFree(_devMIValues))
        throw CUDAError();
    if (cudaSuccess != cudaFree(_devMiIds))
        throw CUDAError();
    if (cudaSuccess != cudaFree(_devIds))
        throw CUDAError();
    if (cudaSuccess != cudaFree(_devDoubleTables))
        throw CUDAError();

    for (int i = 0; i < block_size; i++) {
        _tables[i].finalize();
    }
    delete[] _tables;
    delete[] _hostMIValues;
    delete[] _hostMiIds;

    long myTotalAnal = 0;
    for (auto p : pairs) {
        myTotalAnal += _numSNPs - p.y - 1;
    }
    return myTotalAnal;
}

void MutualInformation::_findNHighestMI(uint3 *_hostMiIds, float *_hostMIValues, uint64_t totalValues, float &minMI,
                                        uint16_t &minMIPos, uint16_t &numEntriesWithMI,
                                        size_t num_outputs, Position *output) {
    int iter = 0;
    Position *auxMI;
    float auxValue;

    if (numEntriesWithMI == 0) { // The first values are directly stored
        for (iter = 0; iter < num_outputs; iter++) {
            auxValue = _hostMIValues[iter];

            auxMI = &output[iter];
            auxMI->rank = auxValue;
            auxMI->p1 = _hostMiIds[iter].x;
            auxMI->p2 = _hostMiIds[iter].y;
            auxMI->p3 = _hostMiIds[iter].z;

            if (auxValue < minMI) {
                minMI = auxValue;
                minMIPos = iter;
            }
        }
        numEntriesWithMI += num_outputs;
    }

    for (; iter < totalValues; iter++) {
        auxValue = _hostMIValues[iter];

        if (auxValue > minMI) { // The value must be inserted
            auxMI = &output[minMIPos];
            auxMI->p1 = _hostMiIds[iter].x;
            auxMI->p2 = _hostMiIds[iter].y;
            auxMI->p3 = _hostMiIds[iter].z;
            auxMI->rank = auxValue;

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
    }
}