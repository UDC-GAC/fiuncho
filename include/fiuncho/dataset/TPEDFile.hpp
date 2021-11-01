/**
 * @file TPEDFile.hpp
 * @date 01/11/2021
 */

#ifndef FIUNCHO_TPEDFILE_HPP
#define FIUNCHO_TPEDFILE_HPP

#include <algorithm>
#include <fiuncho/Dataset.h>
#include <fstream>
#include <string>

namespace TPEDFile
{

/**
 * @brief Class representing a variant from a TPED file, as defined by PLINK
 * (https://www.cog-genomics.org/plink/1.9/formats#tped).
 */

class Variant
{
  public:
    std::string chr_code; // Chromosome code
    std::string v_id;     // Variant identifier
    double pos;           // Position in morgans or centimorgans
    unsigned int coord;   // Base-pair coordinate
    std::vector<char> alleles;

    friend std::istream &operator>>(std::istream &str, Variant &v)
    {
        std::string line;
        std::vector<char> nucleotides;
        Variant tmp;
        if (std::getline(str, line)) {
            std::stringstream iss(line);
            // Parse variant information
            if (!(iss >> tmp.chr_code && iss >> tmp.v_id && iss >> tmp.pos &&
                  iss >> tmp.coord)) {
                throw std::runtime_error("invalid loci information");
            }

            // Read nucleotides
            char c;
            while (iss >> c && (c == 'A' || c == 'C' || c == 'T' || c == 'G')) {
                tmp.alleles.push_back(c);
            }
            if (!iss.eof()) {
                throw std::runtime_error(
                    std::string("invalid nucleotide value '") + c +
                    "' at position " + std::to_string(iss.tellg()));
            }
            if (tmp.alleles.size() % 2) {
                throw std::runtime_error("odd number of nucleotides");
            }
            v.swap(tmp);
        }
        return str;
    }

    void swap(Variant &other)
    {
        std::swap(chr_code, other.chr_code);
        std::swap(v_id, other.v_id);
        std::swap(pos, other.pos);
        std::swap(coord, other.coord);
        std::swap(alleles, other.alleles);
    }
};

/**
 * @brief Class representing a sample from a TFAM file, as defined by PLINK
 * (https://www.cog-genomics.org/plink/1.9/formats#tfam)
 */

class Sample
{
  public:
    std::string fid;   // Family ID
    std::string iid;   // Within-family ID (cannot be '0')
    std::string f_iid; // Within-family ID of father ('0' if father isn't in
                       // the dataset)
    std::string m_iid; // Within-family ID of mother ('0' if mother isn't in
                       // the dataset)
    int sex;           // Sex code ('1' = male, '2' = female, '0' = unknown)
    int ph; // Phenotype value ('1' = control, '2' = case, '-9'/'0'/ //
            // non-numeric = missing data if case/control)

    friend std::istream &operator>>(std::istream &str, Sample &s)
    {
        std::string line;
        Sample tmp;
        if (std::getline(str, line)) {
            std::stringstream iss(line);
            if (iss >> tmp.fid && iss >> tmp.iid && iss >> tmp.f_iid &&
                iss >> tmp.m_iid && iss >> tmp.sex && iss >> tmp.ph &&
                (tmp.ph == 1 || tmp.ph == 2)) {
                /* OK: All read operations worked */
                s.swap(tmp); // C++03 as this answer was written long ago.
            } else {
                throw std::runtime_error("parsing error at " +
                                         std::to_string(iss.tellg()));
            }
        }
        return str;
    }

    void swap(Sample &other)
    {
        std::swap(fid, other.fid);
        std::swap(iid, other.iid);
        std::swap(f_iid, other.f_iid);
        std::swap(m_iid, other.m_iid);
        std::swap(sex, other.sex);
        std::swap(ph, other.ph);
    }
};

void read_samples(const std::string &tfam, std::vector<Sample> &samples,
                  size_t &cases, size_t &ctrls)
{
    std::ifstream file;
    file.open(tfam.c_str(), std::ios::in);
    if (!file.is_open()) {
        throw std::runtime_error("Error while opening " + tfam +
                                 ", check file path/permissions");
    }
    try {
        Sample ind;
        ctrls = 0;
        while (file >> ind) {
            samples.push_back(ind);
            ctrls += ind.ph == 1;
        }
        cases = samples.size() - ctrls;
    } catch (const std::runtime_error &e) {
        throw std::runtime_error("Error in " + tfam + ":" +
                                 std::to_string(samples.size() + 1) + ": " +
                                 e.what());
    }
    file.close();
};

void read_variants(const std::string &tped, const std::vector<Sample> &samples,
                   std::vector<Variant> &variants)
{
    std::ifstream file;
    file.open(tped.c_str(), std::ios::in);
    if (!file.is_open()) {
        throw std::runtime_error("Error while opening " + tped +
                                 ", check file path/permissions");
    }
    try {
        Variant v;
        while (file >> v) {
            if (v.alleles.size() == 2 * samples.size()) {
                variants.push_back(v);
            } else {
                throw std::runtime_error(
                    "Error in " + tped + ":" +
                    std::to_string(variants.size() + 1) +
                    ": the number of nucleotides does not match "
                    "the number of samples");
            }
        }
    } catch (const std::runtime_error &e) {
        throw std::runtime_error("Error in " + tped + ":" +
                                 std::to_string(variants.size() + 1) + ": " +
                                 e.what());
    }
    file.close();
};

template <typename T>
void fill(std::unique_ptr<GenotypeTable<T>[]> &array,
          const std::vector<Sample> &samples,
          const std::vector<Variant> &variants, const size_t cases_words,
          const size_t ctrls_words)
{
    constexpr size_t BITS = sizeof(T) * 8; // Number of bits in T
    // Buffers
    T cases_buff[3], ctrls_buff[3];
    for (size_t i = 0; i < variants.size(); i++) {
        // Find minor allele
        std::map<char, size_t> count;
        std::for_each(variants[i].alleles.begin(), variants[i].alleles.end(),
                      [&count](const char &c) {
                          count.find(c) == count.end() ? count[c] = 1
                                                       : count[c] += 1;
                      });
        char a1 = count.begin()->first;
        char a2 = (++count.begin())->first;
        char minor = count[a1] > count[a2] ? a2 : a1;
        // Clear buffers
        for (auto k = 0; k < 3; k++) {
            cases_buff[k] = 0;
            ctrls_buff[k] = 0;
        }
        // Obtain reference to the next genotype table
        auto &table = array[i];
        // Populate table with the variant information
        size_t cases_cnt = 0;
        size_t ctrls_cnt = 0;
        for (size_t j = 0; j < samples.size(); j++) {
            // Transform alleles into minor allele count
            short c = (variants[i].alleles[j * 2] == minor) +
                      (variants[i].alleles[j * 2 + 1] == minor);
            // For each sample, check phenotype class
            if (samples[j].ph == 1) { // If it's a control append genotype to
                // the 3 control buffers
                for (auto k = 0; k < 3; k++) {
                    ctrls_buff[k] = (ctrls_buff[k] << 1) + (c == k);
                }
                ctrls_cnt++;
            } else { // Else append genotype to the 3 cases buffers
                for (auto k = 0; k < 3; k++) {
                    cases_buff[k] = (cases_buff[k] << 1) + (c == k);
                }
                cases_cnt++;
            }
            // If the buffer is full, write buffer into the bit table and
            // clear the buffer
            if (cases_cnt > 0 && cases_cnt % BITS == 0) {
                const int offset = cases_cnt / BITS - 1;
                for (auto k = 0; k < 3; k++) {
                    table.cases[k * cases_words + offset] = cases_buff[k];
                    cases_buff[k] = 0;
                }
            }
            // Do the same for controls
            if (ctrls_cnt > 0 && ctrls_cnt % BITS == 0) {
                const int offset = ctrls_cnt / BITS - 1;
                for (auto k = 0; k < 3; k++) {
                    table.ctrls[k * ctrls_words + offset] = ctrls_buff[k];
                    ctrls_buff[k] = 0;
                }
            }
        }
        // If the number of controls is not divisible by the bits in T
        if (cases_cnt % BITS != 0) {
            const int offset = cases_cnt / BITS;
            for (auto k = 0; k < 3; k++) {
                table.cases[k * cases_words + offset] = cases_buff[k];
            }
        }
        // Repeat for ctrls
        if (ctrls_cnt % BITS != 0) {
            // Write last (incomplete) word from each row of the table
            const int offset = ctrls_cnt / BITS;
            for (auto k = 0; k < 3; k++) {
                table.ctrls[k * ctrls_words + offset] = ctrls_buff[k];
            }
        }
        // Write 0 in the remaining uninitialized words of the controls
        // array
        for (auto i = (cases_cnt + BITS - 1) / BITS; i < cases_words; i++) {
            for (auto k = 0; k < 3; k++) {
                table.cases[k * cases_words + i] = 0;
            }
        }
        // Repeat for ctrls
        for (auto i = (ctrls_cnt + BITS - 1) / BITS; i < ctrls_words; i++) {
            for (auto k = 0; k < 3; k++) {
                table.ctrls[k * ctrls_words + i] = 0;
            }
        }
    }
};

template <typename T>
void read(const std::string tped, const std::string tfam, size_t &cases,
          size_t &ctrls, size_t &snps,
          std::unique_ptr<GenotypeTable<T>[]> &array)
{
    std::vector<Sample> samples;
    std::vector<Variant> variants;
    // Read input files
    read_samples(tfam, samples, cases, ctrls);
    read_variants(tped, samples, variants);
    snps = variants.size();
#ifdef ALIGN
    // Pad table rows so that each row is divisible by the VPU width
    constexpr size_t NT = ALIGN / sizeof(T); // Number of T's in ALIGN bytes
    constexpr size_t NBITS = ALIGN * 8;      // Number of bits in ALIGN bytes
    const size_t cases_words = (cases + NBITS - 1) / NBITS * NT,
                 ctrls_words = (ctrls + NBITS - 1) / NBITS * NT;
#else
    constexpr size_t NBITS = sizeof(T) * 8; // Number of bits in T
    const size_t cases_words = (cases + NBITS - 1) / NBITS,
                 ctrls_words = (ctrls + NBITS - 1) / NBITS;
#endif
    array = GenotypeTable<T>::make_array(snps, 1, cases_words, ctrls_words);
    // Fill array with the information
    fill(array, samples, variants, cases_words, ctrls_words);
};

} // namespace TPEDFile

#endif
