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
 * @file cpu/CPUEngine.h
 * @author Christian Ponte
 * @date 18 March 2020
 *
 * @brief CPUEngine class declaration, implementing the abstract class Engine.
 */

#ifndef FIUNCHO_CPUENGINE_H
#define FIUNCHO_CPUENGINE_H

#include <fiuncho/dataset/Dataset.h>
#include <fiuncho/engine/Engine.h>
#include <fiuncho/engine/Position.h>

class CPUEngine : public Engine {
  public:
    CPUEngine(int num_proc, int proc_id, int num_threads, bool use_mi,
              Statistics &statistics);

    virtual ~CPUEngine() = default;

    void run(std::string tped, std::string tfam,
             std::vector<Position> &mutual_info, size_t num_outputs) override;

  private:
    class Shared_block {
      public:
        Shared_block(int tid, Dataset &dataset, uint16_t numOutputs,
                     Statistics &statistics)
            : tid(tid), dataset(dataset), numOutputs(numOutputs),
              positions(new Position[numOutputs]), statistics(statistics) {}

        ~Shared_block() { delete[] positions; }

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

#endif
