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

#include <fiuncho/engine/cpu/ContingencyTable.h>

ContingencyTable::ContingencyTable(std::vector<uint32_t> *ctrls1,
                                   std::vector<uint32_t> *cases1,
                                   std::vector<uint32_t> *ctrls2,
                                   std::vector<uint32_t> *cases2)
    : cells(9), ctrl_table(new uint32_t[cells * ctrls1[0].size()]),
      ctrl_words(ctrls1[0].size()),
      case_table(new uint32_t[cells * cases1[0].size()]),
      case_words(cases1[0].size()) {
    for (short i = 0; i < 3; i++) {
        for (short j = 0; j < 3; j++) {
            for (auto word = 0; word < ctrl_words; word++) {
                ctrl_table[(i * 3 + j) * ctrl_words + word] =
                    ctrls1[i][word] & ctrls2[j][word];
            }
        }
    }

    for (short i = 0; i < 3; i++) {
        for (short j = 0; j < 3; j++) {
            for (auto word = 0; word < case_words; word++) {
                case_table[(i * 3 + j) * case_words + word] =
                    cases1[i][word] & cases2[j][word];
            }
        }
    }
}

ContingencyTable::ContingencyTable(const ContingencyTable &minor,
                                   std::vector<uint32_t> *ctrls,
                                   std::vector<uint32_t> *cases)
    : cells(minor.cells * 3), ctrl_words(ctrls[0].size()),
      ctrl_table(new uint32_t[cells * ctrl_words]), case_words(cases[0].size()),
      case_table(new uint32_t[cells * case_words]) {
    for (short i = 0; i < minor.cells; i++) {
        for (short j = 0; j < 3; j++) {
            for (auto word = 0; word < ctrl_words; word++) {
                ctrl_table[(i * 3 + j) * ctrl_words + word] =
                    minor.ctrl_table[i * ctrl_words + word] & ctrls[j][word];
            }
        }
    }

    for (short i = 0; i < minor.cells; i++) {
        for (short j = 0; j < 3; j++) {
            for (auto word = 0; word < case_words; word++) {
                case_table[(i * 3 + j) * case_words + word] =
                    minor.case_table[i * case_words + word] & cases[j][word];
            }
        }
    }
}

ContingencyTable::~ContingencyTable() {
    delete[] ctrl_table;
    delete[] case_table;
}

std::size_t ContingencyTable::get_cells() const { return cells; }

uint16_t ContingencyTable::count_controls(std::size_t cell) const {
    uint16_t count = 0;
    for (auto word = 0; word < ctrl_words; word++) {
        count += popcount(ctrl_table[cell * ctrl_words + word]);
    }
    return count;
}

uint16_t ContingencyTable::count_cases(std::size_t cell) const {
    uint16_t count = 0;
    for (auto word = 0; word < case_words; word++) {
        count += popcount(case_table[cell * case_words + word]);
    }
    return count;
}

uint32_t ContingencyTable::popcount(uint32_t v) const {
    uint32_t u;
    u = v - ((v >> 1) & 033333333333) - ((v >> 2) & 011111111111);
    return ((u + (u >> 3)) & 030707070707) % 63;
}