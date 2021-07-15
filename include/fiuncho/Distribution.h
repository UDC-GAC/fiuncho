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
 * @file Distribution.h
 * @author Christian Ponte
 */

#ifndef FIUNCHO_DISTRIBUTION_H
#define FIUNCHO_DISTRIBUTION_H

#include <type_traits>
#include <vector>

/**
 * @class Distribution
 * @brief Class that represents the distribution of the workload among all the
 * compute units available. The distribution creates an assignment of all
 * *k*-combinations without repetition from a set of *n* SNPs following a
 * Round-robin distribution. To do this, combinations are enumerated using a
 * particular *step* size and an initial *offset*.
 *
 * @tparam T Data type used to index the SNPs in the data set. Needs to be a
 * _signed_ type
 */

template <typename T> class Distribution
{
  public:
    /**
     * @name Attributes
     */
    //@{

    /**
     * Number of SNPs in the set
     */
    const T n;

    /**
     * Size of the combinations
     */
    const T k;

    /**
     * Number of combinations to advance in each iteration
     */
    const T step;

    /**
     * Initial offset when iterating over the combinations
     */
    const T offset;

    //@}

    class const_iterator
    {
        friend class Distribution<T>;
        const T n, k, step, offset;
        std::vector<T> c;

        inline void increment_left()
        {
            // p needs to be a signed value
            T p = k - 2;
            // find first index to increment
            while (p >= 0 && c[p] + 1 > n - (k - p)) {
                --p;
            }
            // if we moved past the first index, that means the list of
            // combinations is past its end
            if (p < 0) {
                c[0] = n;
            } else {
                ++c[p];
                for (auto i = p + 1; i < k - 1; ++i) {
                    c[i] = c[i - 1] + 1;
                }
            }
        }

        inline void increment(T x)
        {
            c[k - 1] += x;
            while (c[0] < n && c[k - 1] >= n) {
                increment_left();
                c[k - 1] += (c[k - 2] + 1) - n;
            }
        }

      public:
        using value_type = T;
        using reference = T;
        using iterator_category = std::input_iterator_tag;
        using pointer = T *;
        using difference_type = void;

        const_iterator(const Distribution<T> &d)
            : n(d.n), k(d.k), step(d.step), offset(d.offset), c(k)
        {
            for (auto i = 0; i < k; ++i) {
                c[i] = i;
            }
            increment(offset);
        }

        const std::vector<T> &operator*() const { return c; }

        const std::vector<T> *operator->() const { return &c; }

        const T &operator[](size_t idx) const { return c[idx]; }

        const_iterator &operator++()
        { // preincrement
            increment(step);
            return *this;
        }

        bool operator<(const const_iterator &val)
        {
            auto i = c.begin();
            auto j = val.c.begin();
            for (; i < c.end() && j < val.c.end(); ++i, ++j) {
                if (*i < *j) {
                    return true;
                } else if (*i > *j) {
                    return false;
                }
            }
            return false;
        }
    };

    /**
     * @name Constructors
     */
    //@{

    /**
     * Create a new distribution.
     *
     * @param n Number of SNPs in the set
     * @param k Size of the combinations to consider
     * @param step Number of combinations to advance between iterations
     * @param offset Number of combinations to advance before starting the
     * iteration
     */

    Distribution(const T &n, const T &k, const T &step, const T &offset)
        : n(n), k(k), step(step), offset(offset)
    {
        static_assert(std::is_signed<T>::value,
                      "Distribution template parameter requires a signed type");
    };

    //@}

    /**
     * @name Methods
     */
    //@{

    /**
     * Create a new distribution from an existing one by layering a secondary
     * *step* and *offset*. The new distribution will conserve the same number
     * of *n* SNPs in the set and combination size *k*. The new *step* is
     * calculated as \f$ step_{1} * step_{2} \f$, and the new *offset* is
     * calculated as \f$ offset_{1} * step_{2} + offset_{2} \f$.
     *
     * @param step Secondary step to include in the new distribution
     * @param offset Secondary offset to include in the new distribution
     */

    Distribution<T> layer(const T step, const T offset) const
    {
        return Distribution<T>(n, k, this->step * step,
                               this->offset * step + offset);
    }

    /**
     * Returns an iterator to the first combination of the distribution.
     * Combinations returned reuse the same std::vector object, succesively
     * updating its content.
     *
     * @return Iterator to the first combination
     */

    const_iterator begin() const { return const_iterator(*this); }

    /**
     * Returns an iterator to the combination following the last combination of
     * the distribution.
     *
     * @return Iterator to the combination following the last combination.
     */

    const_iterator end() const
    {
        const_iterator end_iter(*this);
        end_iter.c[0] = n;
        return end_iter;
    }

    //@}
};

#endif
