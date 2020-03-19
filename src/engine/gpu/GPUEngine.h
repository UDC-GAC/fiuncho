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
 * @file gpu/GPUEngine.h
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief GPUEngine class declaration.
 */

#ifndef MPI3SNP_GPUENGINE_H
#define MPI3SNP_GPUENGINE_H

#include "Engine.h"
#include "Dataset.h"
#include "Position.h"

class GPUEngine : public Engine {
public:
    GPUEngine(unsigned int proc_num, unsigned int proc_id, std::vector<std::pair<unsigned int, unsigned int>> gpu_map,
              bool use_mi, Statistics &statistics);

    void run(std::string tped, std::string tfam, std::vector<Position> &mutual_info, size_t num_outputs) override;

private:
    unsigned int proc_num, proc_id, gpu_id;
    bool use_mi;
    Statistics &statistics;
};

#endif //MPI3SNP_GPUENGINE_H
