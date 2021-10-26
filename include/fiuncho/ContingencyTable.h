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

    ContingencyTable() : size(0), alloc(nullptr), cases(nullptr), ctrls(nullptr)
    {
    }

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
          alloc(new T[size * 2 + N], std::default_delete<T[]>()),
#undef N
#undef NBITS
          cases((T *)((((uintptr_t)alloc.get()) + ALIGN - 1) / ALIGN * ALIGN)),
#else
          size(std::pow(3, order)),
          alloc(new T[size * 2], std::default_delete<T[]>()),
          cases(alloc.get()),
#endif
          ctrls(cases + size){};

    static std::unique_ptr<ContingencyTable<T>[]>
    make_array(const size_t N, const short order) {
    // Allocate a single contiguous array for all subtables
#ifdef ALIGN
        constexpr size_t NT = ALIGN / sizeof(T);
        const size_t size = (size_t)(std::pow(3, order) + NT - 1) / NT * NT;
        std::shared_ptr<T> alloc(new T[N * size * 2 + NT],
                                 std::default_delete<T[]>());
        T *ptr = (T *)((((uintptr_t)alloc.get()) + ALIGN - 1) / ALIGN * ALIGN);
#else
        const size_t size = std::pow(3, order);
        std::shared_ptr<T> alloc(new T[N * size * 2],
                                 std::default_delete<T[]>());
        T *ptr = alloc.get();
#endif
        // Allocate the ContingencyTable object array
        auto array = std::make_unique<ContingencyTable<T>[]>(N);
        // Initialize ContingencyTable objects
        for (size_t i = 0; i < N; ++i) {
            new (array.get() + i) ContingencyTable(
                size, alloc, ptr + i * size * 2, ptr + i * size * 2 + size);
        }
        return array;
    }

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
    ContingencyTable(const size_t size, std::shared_ptr<T> alloc, T *cases,
                     T *ctrls)
        : size(size), alloc(alloc), cases(cases), ctrls(ctrls)
    {
    }

    std::shared_ptr<T> alloc;

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
