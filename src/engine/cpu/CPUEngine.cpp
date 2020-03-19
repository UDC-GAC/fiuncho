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
 * @file cpu/CPUEngine.cpp
 * @author Christian Ponte
 * @date 18 March 2020
 *
 * @brief CPUEngine class members implementation.
 */

#include "CPUEngine.h"
#include "MutualInformation.h"
#include "ThreadError.h"
#include "Distributor.h"
#include <thread>

CPUEngine::CPUEngine(int num_proc, int proc_id, int num_threads, bool use_mi, Statistics &statistics) :
        num_proc(num_proc),
        proc_id(proc_id),
        num_threads(num_threads),
        use_mi(use_mi),
        statistics(statistics) {}

void CPUEngine::run(std::string tped, std::string tfam, std::vector<Position> &mutual_info, size_t num_outputs) {
    statistics.Begin_timer("SNPs read time");
    Dataset *dataset;
    try {
        dataset = new Dataset(tped, tfam, Dataset::Regular);
    } catch (const Dataset::Read_error &error) {
        throw Engine::Error(error.what());
    }
    statistics.End_timer("SNPs read time");

    Distributor<uint32_t, std::pair<uint32_t, uint32_t>> distributor(dataset->get_SNP_count(), num_proc * num_threads);
    statistics.Addi("SNP count", dataset->get_SNP_count());
    statistics.Addi("Number of cases", dataset->get_case_count());
    statistics.Addi("Number of controls", dataset->get_ctrl_count());

    std::vector<std::thread *> threads;
    std::vector<Shared_block *> params;
    
    // Create thread entities that call to the functions below
    for (int tid = 0; tid < num_threads; tid++) {
        params.push_back(new Shared_block( tid, *dataset, num_outputs, statistics));
        distributor.get_pairs(
            [](uint32_t first, uint32_t second) {
                return std::make_pair(first, second); 
            },
            proc_id * num_threads + tid,
            params.back()->pairs
        );
        threads.push_back(new std::thread(CPUEngine::thread, params.back()));
    }

    std::vector<Position> auxMutualInfo(num_threads * num_outputs);

    // Wait for the completion of all threads
    for (auto i = 0; i < threads.size(); i++) {
        threads[i]->join();
        delete threads[i];
        memcpy(&auxMutualInfo[i * num_outputs], params[i]->positions, num_outputs * sizeof(Position));
        delete params[i];
    }

    delete dataset;

    // Sort the auxiliar array and print the results
    std::sort(auxMutualInfo.begin(), auxMutualInfo.end());
    mutual_info.resize(num_outputs);
    memcpy(&mutual_info[0], &auxMutualInfo[num_outputs * (num_threads - 1)], sizeof(Position) * num_outputs);
}

void CPUEngine::thread(CPUEngine::Shared_block *params) {
    Algorithm<std::pair<uint32_t, uint32_t>> *search = new MutualInformation(
            params->dataset.get_SNP_count(), params->dataset.get_case_count(), params->dataset.get_cases(),
            params->dataset.get_ctrl_count(), params->dataset.get_ctrls());

    // Variables to work with the outputs
    Position *positions = new Position[params->numOutputs];

    std::string timer_label;
    timer_label += "Thread " + std::to_string(params->tid) + " runtime";
    std::string analysis_label;
    analysis_label += "Thread " + std::to_string(params->tid) + " computations";

    params->statistics.Begin_timer(timer_label);

    long myTotalAnal = search->compute(params->pairs, params->numOutputs, positions);

    params->statistics.End_timer(timer_label);
    params->statistics.Addl(analysis_label, myTotalAnal);

    memcpy(params->positions, positions, params->numOutputs * sizeof(Position));

    delete[] positions;
    delete search;
}
