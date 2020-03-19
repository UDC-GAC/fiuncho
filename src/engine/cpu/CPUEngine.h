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
 * @file cpu/CPUEngine.h
 * @author Christian Ponte
 * @date 18 March 2020
 *
 * @brief CPUEngine class declaration, implementing the abstract class Engine.
 */

#ifndef MPI3SNP_CPUENGINE_H
#define MPI3SNP_CPUENGINE_H

#include "Engine.h"
#include "Dataset.h"
#include "Position.h"

class CPUEngine : public Engine {
public:
    CPUEngine(int num_proc, int proc_id, int num_threads, bool use_mi, Statistics &statistics);

    virtual ~CPUEngine() = default;

    void run(std::string tped, std::string tfam, std::vector<Position> &mutual_info, size_t num_outputs) override;

private:
    class Shared_block {
    public:
        Shared_block(int tid, Dataset &dataset, uint16_t numOutputs, Statistics &statistics) :
                tid(tid),
                dataset(dataset),
                numOutputs(numOutputs),
                positions(new Position[numOutputs]),
                statistics(statistics) {}

        ~Shared_block() {
            delete[] positions;
        }

        const int tid;
        Dataset &dataset;
        std::vector<std::pair<uint32_t, uint32_t>> pairs;
        const uint16_t numOutputs;
        Position *positions;
        Statistics &statistics;
    };

    static void thread(CPUEngine::Shared_block *params);

    int num_proc;
    int proc_id;
    int num_threads;
    bool use_mi;
    Statistics &statistics;
};

#endif //MPI3SNP_CPUENGINE_H
