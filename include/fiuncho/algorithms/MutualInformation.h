/**
 * @file MutualInformation.h
 * @author Christian Ponte
 */

#ifndef FIUNCHO_MUTUALINFORMATION_H
#define FIUNCHO_MUTUALINFORMATION_H

#include <fiuncho/algorithms/Algorithm.h>

#include <cmath>

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

template <class T> class MutualInformation : public Algorithm<T>
{
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

    MutualInformation(unsigned int num_cases, unsigned int num_ctrls)
    {
        inv_inds = 1.0 / (num_cases + num_ctrls);
        float p = num_cases * inv_inds;
        h_y = (-1.0) * p * logf(p);
        p = num_ctrls * inv_inds;
        h_y -= p * logf(p);
    };

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
