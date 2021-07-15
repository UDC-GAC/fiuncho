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
 * @file MutualInfo.h
 * @author Christian Ponte
 * @date 30 May 2018
 *
 * @brief MutualInfo structure definition, used for the mutual information
 * outputs with the result and the SNPs associated.
 */

#ifndef FIUNCHO_RESULT_H
#define FIUNCHO_RESULT_H

#include <cstring>
#include <sstream>
#include <vector>

template <typename U, typename V> class Result
{
  public:
    Result() : val(0) {}

    bool operator<(const Result &rhs) const { return val < rhs.val; }

    bool operator>(const Result &rhs) const { return rhs < *this; }

    bool operator<=(const Result &rhs) const { return !(rhs < *this); }

    bool operator>=(const Result &rhs) const { return !(*this < rhs); }

    bool operator==(const Result &rhs) const
    {
        return val == rhs.val && combination == rhs.combination;
    }

    bool operator!=(const Result &rhs) const { return !(*this == rhs); }

    inline std::string str() const
    {
        std::stringstream buffer;
        for (U i : combination) {
            buffer << i << " ";
        }
        buffer << val;
        return buffer.str();
    }

    static std::ostream &serialize(std::ostream &os, const Result<U, V> &r)
    {
        auto size = r.combination.size();
        os.write(reinterpret_cast<char const *>(&size), sizeof(size));
        os.write(reinterpret_cast<char const *>(r.combination.data()),
                 size * sizeof(U));
        os.write(reinterpret_cast<char const *>(&r.val), sizeof(V));
        return os;
    }

    static std::istream &deserialize(std::istream &is, Result &r)
    {
        decltype(r.combination.size()) size;
        is.read(reinterpret_cast<char *>(&size), sizeof(size));
        r.combination.resize(size);
        is.read(reinterpret_cast<char *>(r.combination.data()),
                size * sizeof(U));
        is.read(reinterpret_cast<char *>(&r.val), sizeof(V));
        return is;
    }

    std::vector<U> combination;
    V val;
};

#endif
