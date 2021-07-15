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

#include <cstdint>
#include <fiuncho/ThreadedSearch.h>
#include <fiuncho/dataset/Dataset.h>
#include <fiuncho/utils/MaxArray.h>
#include <iostream>
#include <thread>
#include <time.h>
#include <vector>

/*
 *  This benchmark programs calculates bit tables of a specific order following
 *  this process:
 *      1. Spawn as many threads as indicated
 *      2. Set processor affinity as indicated
 *      3. Initialize previous tables and such, and synchronize threads at the
 *         end of the initialization
 *      4. Warm up the CPU core by running an additional 10% of the total
 *         iterations
 *      5. Measure bit table computation time
 *      6. Print elapsed time on each thread
 *
 *  Program arguments:
 *      1: Comma-separated list of hardware threads to run on.
 *      2: Order of the tables to benchmark
 *      3: How many times the main loop is repeated
 *      4: Path to the TPED input file
 *      5: Path to the TFAM input file
 */

int main(int argc, char *argv[])
{
    if (argc != 5) {
        std::cout << argv[0] << " <NTHREADS> <ORDER> <TPED> <TFAM>"
                  << std::endl;
        return 0;
    }

    // Initialization
    // Arguments
    const unsigned short thread_count = atoi(argv[1]);
    const unsigned short order = atoi(argv[2]);
    const std::string tped = argv[3], tfam = argv[4];
    // Data
#ifdef ALIGN
    const auto dataset = Dataset<uint64_t>::read<ALIGN>(tped, tfam);
#else
    const auto dataset = Dataset<uint64_t>::read(tped, tfam);
#endif
    Distribution<int> distribution(dataset.snps, order - 1, 1, 0);
    auto search = ThreadedSearch(thread_count);
    search.run(dataset, order, distribution, 10);

    return 0;
}
