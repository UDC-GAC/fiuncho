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
 * @file ContingencyTable.cpp
 * @author Christian Ponte
 * @date 29 May 2018
 *
 * @brief ContingencyTable class members implementation.
 */

#include <cmath>
#include <fiuncho/ContingencyTable.h>

template <>
ContingencyTable<uint32_t>::ContingencyTable(const short order,
                                             const size_t cases_words,
                                             const size_t ctrls_words)
    : size(std::pow(3, order)), cases_words(cases_words),
      ctrls_words(ctrls_words), alloc(std::make_unique<uint32_t[]>(size * 2)),
      cases(alloc.get()), ctrls(cases + size)
{
}
