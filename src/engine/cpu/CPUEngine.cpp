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
 * @file cpu/CPUEngine.cpp
 * @author Christian Ponte
 * @date 18 March 2020
 *
 * @brief CPUEngine class members implementation.
 */

#include <fiuncho/engine/Distributor.h>
#include <fiuncho/engine/cpu/CPUEngine.h>
#include <fiuncho/engine/cpu/MutualInformation.h>
#include <fiuncho/engine/cpu/ThreadError.h>
#include <thread>

CPUEngine::CPUEngine(int num_proc, int proc_id, int num_threads, bool use_mi,
                     Statistics &statistics)
    : num_proc(num_proc), proc_id(proc_id), num_threads(num_threads),
      use_mi(use_mi), statistics(statistics) {}

void CPUEngine::run(std::string tped, std::string tfam,
                    std::vector<Position> &mutual_info, size_t num_outputs) {
    statistics.begin_timer("SNPs read time");
    Dataset *dataset;
    try {
        dataset = new Dataset(tped, tfam, Dataset::Regular);
    } catch (const Dataset::Read_error &error) {
        throw Engine::Error(error.what());
    }
    statistics.end_timer("SNPs read time");

    Distributor<uint32_t, std::pair<uint32_t, uint32_t>> distributor(
        dataset->get_SNP_count(), num_proc * num_threads);
    statistics.addi("SNP count", dataset->get_SNP_count());
    statistics.addi("Number of cases", dataset->get_case_count());
    statistics.addi("Number of controls", dataset->get_ctrl_count());

    std::vector<std::thread *> threads;
    std::vector<Shared_block *> params;

    // Create thread entities that call to the functions below
    for (int tid = 0; tid < num_threads; tid++) {
        params.push_back(
            new Shared_block(tid, *dataset, num_outputs, statistics));
        distributor.get_pairs(
            [](uint32_t first, uint32_t second) {
                return std::make_pair(first, second);
            },
            proc_id * num_threads + tid, params.back()->pairs);
        threads.push_back(new std::thread(CPUEngine::thread, params.back()));
    }

    std::vector<Position> auxMutualInfo(num_threads * num_outputs);

    // Wait for the completion of all threads
    for (auto i = 0; i < threads.size(); i++) {
        threads[i]->join();
        delete threads[i];
        memcpy(&auxMutualInfo[i * num_outputs], params[i]->positions,
               num_outputs * sizeof(Position));
        delete params[i];
    }

    delete dataset;

    // Sort the auxiliar array and print the results
    std::sort(auxMutualInfo.begin(), auxMutualInfo.end());
    mutual_info.resize(num_outputs);
    memcpy(&mutual_info[0], &auxMutualInfo[num_outputs * (num_threads - 1)],
           sizeof(Position) * num_outputs);
}

void CPUEngine::thread(CPUEngine::Shared_block *params) {
    Algorithm<std::pair<uint32_t, uint32_t>> *search = new MutualInformation(
        params->dataset.get_SNP_count(), params->dataset.get_case_count(),
        params->dataset.get_cases(), params->dataset.get_ctrl_count(),
        params->dataset.get_ctrls());

    // Variables to work with the outputs
    Position *positions = new Position[params->numOutputs];

    std::string timer_label;
    timer_label += "Thread " + std::to_string(params->tid) + " runtime";
    std::string analysis_label;
    analysis_label += "Thread " + std::to_string(params->tid) + " computations";

    params->statistics.begin_timer(timer_label);

    long myTotalAnal =
        search->compute(params->pairs, params->numOutputs, positions);

    params->statistics.end_timer(timer_label);
    params->statistics.addl(analysis_label, myTotalAnal);

    memcpy(params->positions, positions, params->numOutputs * sizeof(Position));

    delete[] positions;
    delete search;
}
