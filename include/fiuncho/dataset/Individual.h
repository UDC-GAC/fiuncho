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
 * @file Individual.h
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief Individual class definition and implementation.
 */

#ifndef FIUNCHO_INDIVIDUAL_H
#define FIUNCHO_INDIVIDUAL_H

#include <sstream>
#include <stdexcept>
#include <string>

struct Individual {
    class InvalidIndividual : public std::runtime_error {
      public:
        InvalidIndividual(const std::string &msg) : std::runtime_error(msg){};

        virtual ~InvalidIndividual(){};
    };

    std::string fid;   // Family ID
    std::string iid;   // Within-family ID (cannot be '0')
    std::string f_iid; // Within-family ID of father ('0' if father isn't in
                       // the dataset)
    std::string m_iid; // Within-family ID of mother ('0' if mother isn't in
                       // the dataset)
    int sex;           // Sex code ('1' = male, '2' = female, '0' = unknown)
    int ph; // Phenotype value ('1' = control, '2' = case, '-9'/'0'/ //
            // non-numeric = missing data if case/control)

    friend std::istream &operator>>(std::istream &str, Individual &ind) {
        std::string line;
        Individual tmp;
        if (std::getline(str, line)) {
            std::stringstream iss(line);
            if (iss >> tmp.fid && iss >> tmp.iid && iss >> tmp.f_iid &&
                iss >> tmp.m_iid && iss >> tmp.sex && iss >> tmp.ph &&
                (tmp.ph == 1 || tmp.ph == 2)) {
                /* OK: All read operations worked */
                ind.swap(tmp); // C++03 as this answer was written long ago.
            } else {
                throw InvalidIndividual("parsing error at " +
                                        std::to_string(iss.tellg()));
            }
        }
        return str;
    }

    void swap(Individual &other) {
        std::swap(fid, other.fid);
        std::swap(iid, other.iid);
        std::swap(f_iid, other.f_iid);
        std::swap(m_iid, other.m_iid);
        std::swap(sex, other.sex);
        std::swap(ph, other.ph);
    }
};

#endif
