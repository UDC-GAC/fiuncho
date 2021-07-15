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
 * @file MutualInformation.h
 * @author Christian Ponte
 */

#ifndef FIUNCHO_MUTUALINFORMATION_H
#define FIUNCHO_MUTUALINFORMATION_H

#include <fiuncho/algorithms/Algorithm.h>

/**
 * @class MutualInformation
 * @brief Class implementing the Mutual Information (MI) computation. Taking the
 * genotype and phenotype variability as two random variables \f$X\f$ and
 * \f$Y\f$, MI can be obtained as:
 *
 * \f{equation}{
 * MI(X; Y) = H(X) + H(Y) - H(X,Y)
 * \f}
 *
 * Marginal entropies of a single and two variables are defined as:
 *
 * \f{equation}{
 * H(X) = - \sum_{x \in X} p(x)\log p(x)
 * \f}
 * \f{equation}{
 * H(X,Y) = - \sum_{x, y} p(x, y)\log p(x, y)
 * \f}
 *
 * These probabilities can be obtained directly from a ContingencyTable as
 * the division between the number of occurrences and the number of total
 * observations.
 *
 * @tparam T data type used to represent the MI values
 */

template <class T>
class MutualInformation : public Algorithm<T> {
  public:
    /**
     * @name Constructors
     */
    //@{

    /**
     * Create a MutualInformation instance that allows to compute the MI value
     * of a particular ContingencyTable for a fixed number of cases and
     * controls. The number of cases and controls is fixed so that multiple
     * calculations of the MI for the same data set do not require to compute
     * \f$H(Y)\f$ repeatedly.
     *
     * @param num_cases Path to the tped data file
     * @param num_ctrls Path to the tfam data file
     */

    MutualInformation(unsigned int num_cases, unsigned int num_ctrls);

    //@}

    /**
     * @name Methods
     */
    //@{

    /**
     * Compute the MI of a particular ContingencyTable.
     *
     * @return The MI value
     * @param table ContingencyTable from which the MI value is to be calculated
     * @tparam U Data type used in the input ContingencyTable's to represent the
     * count of individuals
     */

    template <class U>
    T compute(const ContingencyTable<U> &table) const noexcept;

    //@}

  private:
    T inv_inds;
    // Entropy of Y
    T h_y;
};

#endif
