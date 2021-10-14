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
 * @file ContingencyTable.h
 * @author Christian Ponte
 */

#ifndef FIUNCHO_CONTINGENCYTABLE_H
#define FIUNCHO_CONTINGENCYTABLE_H

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>

/**
 * @class ContingencyTable
 * @brief Class representing the genotypes frequencies of the combination of \a
 * S SNPs for all individuals. The size of the table is: \f$ size = 3^S \f$.
 *
 * The genotype frequencies are segregated by the class of the individuals
 * (cases or controls) in two subtables.
 *
 * @tparam T data type used to represent the frequencies in the table
 */

template <class T> class ContingencyTable
{
  public:
    ContingencyTable(const ContingencyTable<T> &) = delete;
    ContingencyTable(ContingencyTable<T> &&) = default;

    /**
     * @name Constructors
     */
    //@{

    /**
     * Create a new uninitialized table, allocating an array with enough space
     * to hold \f$ 3^{order} \f$ values of type \a T.
     *
     * @param order Number of SNPs represented in combination inside the table
     */

    ContingencyTable(const short order)
        :
#ifdef ALIGN
#define N (ALIGN / sizeof(T))
#define NBITS (ALIGN * 8)
          size((size_t)(std::pow(3, order) + N - 1) / N * N),
          alloc(std::make_unique<T[]>(size * 2 + N)),
#undef N
#undef NBITS
          cases((T *)((((uintptr_t)alloc.get()) + ALIGN - 1) / ALIGN * ALIGN)),

#else
          size(std::pow(3, order)), alloc(std::make_unique<T[]>(size * 2)),
          cases(alloc.get()),
#endif
          ctrls(cases + size)

              {};

    //@}

    /**
     * @name Attributes
     */
    //@{

    /**
     * Number of values in each subtable
     */
    const size_t size;

  private:
    std::unique_ptr<T[]> alloc;

  public:
    /**
     * Subtable containing the genotype frequencies for the case group
     */
    T *cases;

    /**
     * Subtable containing the genotype frequencies for the control group
     */
    T *ctrls;

    //@}
};

#endif
