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
 * @file SNP.h
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief SNP class definition and implementation.
 */

#ifndef FIUNCHO_SNP_H
#define FIUNCHO_SNP_H

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct SNP {
    class InvalidSNP : public std::runtime_error {
      public:
        InvalidSNP(const std::string &message) : std::runtime_error(message){};

        virtual ~InvalidSNP(){};
    };

    std::string chr_code; // Chromosome code
    std::string v_id;     // Variant identifier
    double pos;           // Position in morgans or centimorgans
    unsigned int coord;   // Base-pair coordinate
    std::vector<uint8_t> genotypes;

    friend std::istream &operator>>(std::istream &str, SNP &snp) {
        std::string line;
        std::vector<char> nucleotides;
        SNP tmp;
        if (std::getline(str, line)) {
            std::stringstream iss(line);
            // Parse SNP information
            if (!(iss >> tmp.chr_code && iss >> tmp.v_id && iss >> tmp.pos &&
                  iss >> tmp.coord)) {
                throw InvalidSNP("invalid loci information");
            }

            // Read nucleotides
            char c;
            while (iss >> c && (c == 'A' || c == 'C' || c == 'T' || c == 'G')) {
                nucleotides.push_back(c);
            }

            if (!iss.eof()) {
                throw InvalidSNP(std::string("invalid nucleotide value '") + c +
                                 "' at position " +
                                 std::to_string(iss.tellg()));
            }

            if (nucleotides.size() % 2) {
                throw InvalidSNP("odd number of nucleotides");
            }

            // Find minor character
            char a1 = nucleotides[0];
            auto pos = std::find_if(nucleotides.begin(), nucleotides.end(),
                                    [&a1](char n) { return n != a1; });
            char a2 = pos == nucleotides.end() ? a1 : *pos;
            a1 = std::min(a1, a2);

            // Translate nucleotides into genotype
            tmp.genotypes.reserve(nucleotides.size() / 2);
            uint8_t allele;
            for (auto it = nucleotides.begin(); it < nucleotides.end();) {
                allele = (*(it++) != a1) + (*(it++) != a1);
                tmp.genotypes.push_back(allele);
            }

            snp.swap(tmp);
        }
        return str;
    }

    void swap(SNP &other) {
        std::swap(chr_code, other.chr_code);
        std::swap(v_id, other.v_id);
        std::swap(pos, other.pos);
        std::swap(coord, other.coord);
        std::swap(genotypes, other.genotypes);
    }
};

#endif
