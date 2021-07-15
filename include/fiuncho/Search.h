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
 * @file Search.h
 * @author Christian Ponte
 */

#include <fiuncho/Distribution.h>
#include <fiuncho/dataset/Dataset.h>
#include <fiuncho/utils/Result.h>

#ifndef FIUNCHO_SEARCH_H
#define FIUNCHO_SEARCH_H

/**
 * Epistasis search class interface. Classes implementing the Search interface
 * make use of the hardware resources available to a particular MPI process to
 * examine the fraction of the SNP combinations assigned to that process.
 */

class Search
{
  public:
    virtual ~Search() {}

    /**
     * @name Methods
     */
    //@{

    /**
     * Run the epistasis search for this MPI process
     *
     * @return Vector of Result's
     * @param dataset Dataset from which the SNPs are read
     * @param order Size of the combinations to explore
     * @param distribution Distribution object for the particular combination
     * distribution of this MPI process
     * @param outputs Number of results to include in the output vector
     */

    virtual std::vector<Result<int, float>>
    run(const Dataset<uint64_t> &dataset, const unsigned short order,
        const Distribution<int> &distribution,
        const unsigned int outputs) = 0;

    //@}
};

#endif
