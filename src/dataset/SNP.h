/*
 * MPI3SNP ~ https://github.com/chponte/mpi3snp
 *
 * Copyright 2018 Christian Ponte
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
 * documentation files (the "Software"), to deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the
 * Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
 * WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * @file SNP.h
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief SNP class definition and implementation.
 */

#ifndef MPI3SNP_SNP_H
#define MPI3SNP_SNP_H

#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

struct SNP {
    class InvalidSNP : public std::runtime_error {
    public:
        InvalidSNP(const std::string &message) : std::runtime_error(message) {};

        virtual ~InvalidSNP() {};
    };

    std::string chr_code; // Chromosome code
    std::string v_id; // Variant identifier
    double pos; // Position in morgans or centimorgans
    unsigned int coord; // Base-pair coordinate
    std::vector<uint8_t> genotypes;

    friend std::istream &operator>>(std::istream &str, SNP &snp) {
        std::string line;
        std::vector<char> nucleotides;
        SNP tmp;
        if (std::getline(str, line)) {
            std::stringstream iss(line);
            // Parse SNP information
            if (!(iss >> tmp.chr_code &&
                  iss >> tmp.v_id &&
                  iss >> tmp.pos &&
                  iss >> tmp.coord)) {
                throw InvalidSNP("invalid loci information");
            }

            // Read nucleotides
            char c;
            while (iss >> c && (c == 'A' || c == 'C' || c == 'T' || c == 'G')) {
                nucleotides.push_back(c);
            }

            if (!iss.eof()) {
                throw InvalidSNP(std::string("invalid nucleotide value '") + c + "' at position " +
                                 std::to_string(iss.tellg()));
            }

            if (nucleotides.size() % 2) {
                throw InvalidSNP("odd number of nucleotides");
            }

            // Find minor character
            char a1 = nucleotides[0];
            auto pos = std::find_if(nucleotides.begin(), nucleotides.end(), [&a1](char n) { return n != a1; });
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

#endif //MPI3SNP_SNP_H
