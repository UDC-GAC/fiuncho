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
 * @file GenotypeTable.h
 * @author Christian Ponte
 */

#ifndef FIUNCHO_GENOTYPETABLE_H
#define FIUNCHO_GENOTYPETABLE_H

#include <fiuncho/ContingencyTable.h>

#include <cstddef>
#include <cstdint>
#include <memory>

/**
 * @class GenotypeTable
 * @brief Class representing the genotypes of \a M individuals using a binary
 * encoding. The size of the table is determined by the number of SNPs \a S
 * represented in combination: \f$ size = 3^S \f$.
 *
 * The genotype information of the individuals is stored in two subtables, one
 * for the individuals from the cases group, and another for the individuals
 * from the controls group.
 *
 * @tparam T data type used to represent the individuals information in the
 * table
 */

template <class T> class GenotypeTable
{
  public:
    GenotypeTable(const GenotypeTable<T> &) = delete;
    GenotypeTable(GenotypeTable<T> &&) = default;

    /**
     * @name Constructors
     */
    //@{

    /**
     * Create a table representing a single SNP, using the array allocations
     * provided.
     *
     * @param cases Array allocation for all rows of the individuals from the
     * cases group
     * @param cases_words Number of values of type \a T required to represent
     * the genotypes for all individuals in the case group
     * @param ctrls Array allocation for all rows of the individuals from the
     * controls group
     * @param ctrls_words Number of values of type \a T required to represent
     * the genotypes for all individuals in the control group
     */

    GenotypeTable(T *cases, const size_t cases_words, T *ctrls,
                  const size_t ctrls_words)
        : order(1), size(3), cases_words(cases_words), ctrls_words(ctrls_words),
          alloc(nullptr), cases(cases), ctrls(ctrls){};

    /**
     * Create a new uninitialized table, allocating an array with enough space
     * to hold \f$ 3^{order} (cases\_words + ctrls\_words) \f$ values of type
     * \a T.
     *
     * @param order Number of SNPs represented in combination inside the table
     * @param cases_words Number of values of type \a T required to represent
     * the genotypes for all individuals in the case group
     * @param ctrls_words Number of values of type \a T required to represent
     * the genotypes for all individuals in the control group
     */

    GenotypeTable(const short order, const size_t cases_words,
                  const size_t ctrls_words);

    //@}

    /**
     * @name Methods
     */
    //@{

    /**
     * Combine the SNPs represented in tables \a t1 and \a t2 into a single
     * table \a out.
     *
     * @param t1 Reference to a GenotypeTable representing any number of SNPs in
     * combination
     * @param t2 Reference to a GenotypeTable representing a single SNP
     * @param out Reference to a GenotypeTable where the combination of \a t1
     * and \a t2 will be stored
     */

    static void combine(const GenotypeTable<T> &t1, const GenotypeTable<T> &t2,
                        GenotypeTable<T> &out) noexcept;

    /**
     * Combine the SNPs represented in tables \a t1 and \a t2, and store the
     * genotype frequencies in the ContingencyTable \a out.
     *
     * @param t1 Reference to a GenotypeTable representing any number of SNPs in
     * combination
     * @param t2 Reference to a GenotypeTable representing a single SNP
     * @param out Reference to a ContingencyTable where the genotype frequencies
     * of the combination of \a t1 and \a t2 will be stored
     */

    template <class U>
    static void combine_and_popcnt(const GenotypeTable<T> &t1,
                                   const GenotypeTable<T> &t2,
                                   ContingencyTable<U> &out) noexcept;

    //@}

    /**
     * @name Attributes
     */
    //@{

    /**
     * Number of SNPs represented in the table
     */
    const size_t order;

    /**
     * Number of rows contained in each subtable
     */
    const size_t size;

    /**
     * Number of values contained in each row of the cases subtable
     */
    const size_t cases_words;

    /**
     * Number of values contained in each row of the controls subtable
     */
    const size_t ctrls_words;

  private:
    std::unique_ptr<T[]> alloc;

  public:
    /**
     * Array representing the cases subtable
     */
    T *cases;

    /**
     * Array representing the controls subtable
     */
    T *ctrls;

    //@}
};

#endif