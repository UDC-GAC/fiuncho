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
 * @date 29 May 2018
 *
 * @brief ContingencyTable class declaration.
 */

#ifndef FIUNCHO_CONTINGENCYTABLE_H
#define FIUNCHO_CONTINGENCYTABLE_H

#include <cstddef>
#include <cstdint>
#include <vector>

class ContingencyTable {
  public:
    ContingencyTable(std::vector<uint32_t> *ctrls1,
                     std::vector<uint32_t> *cases1,
                     std::vector<uint32_t> *ctrls2,
                     std::vector<uint32_t> *cases2);

    ContingencyTable(const ContingencyTable &minor,
                     std::vector<uint32_t> *ctrls,
                     std::vector<uint32_t> *cases);

    ~ContingencyTable();

    std::size_t get_cells() const;

    uint16_t count_controls(std::size_t cell) const;

    uint16_t count_cases(std::size_t cell) const;

  private:
    uint32_t popcount(uint32_t v) const;

    std::size_t cells, ctrl_words, case_words;
    uint32_t *ctrl_table, *case_table;
};

#endif
