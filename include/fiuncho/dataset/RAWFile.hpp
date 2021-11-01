/**
 * @file RAWFile.hpp
 * @date 01/11/2021
 */

#ifndef FIUNCHO_RAWFILE_HPP
#define FIUNCHO_RAWFILE_HPP

#include <fiuncho/Dataset.h>
#include <fiuncho/GenotypeTable.h>
#include <fstream>
#include <regex>
#include <string>

namespace RawFile
{

/**
 * @brief Class representing a single sample from a raw file, implementing the
 * operator>> function for convenient reading. Supports both conventional raw
 * files, as defined by PLINK
 * (https://www.cog-genomics.org/plink/1.9/formats#raw), and GAMETES "raw"
 * files. Relies heavily on the header fields to identify columns to read, and
 * whether or not it is a conventional raw file or if it comes from GAMETES:
 *      1. A sample is from GAMETES if it contains no sample information,
 *      2. Genotype fields contain either no '_', or the '_' is followed by a
 *         {A,C,G,T}.
 *      3. Phenotype column is named "PHENOTYPE" for raw files, and "Class" for
 *         GAMETES.
 *      4. Phenotype for raw files uses '1' for controls, and '2' for cases.
 *         GAMETES uses '0' for controls, and '1' for cases.
 */

class Sample
{
  public:
    Sample() {}

    Sample(bool has_info, const std::vector<bool> &read_mask,
           size_t phenotype_index)
        : has_info(has_info), read_mask(read_mask),
          phenotype_index(phenotype_index)
    {
    }

    /**
     * Family ID
     */
    std::string fid;

    /**
     * Within-family ID (cannot be '0')
     */
    std::string iid;

    /**
     * Within-family ID of father ('0' if father isn't in the dataset)
     */
    std::string f_iid;

    /**
     * Within-family ID of mother ('0' if mother isn't in the dataset)
     */
    std::string m_iid;

    /**
     * Sex code ('1' = male, '2' = female, '0' = unknown)
     */
    short sex;

    /**
     * Phenotype value ('1' = control, '2' = case)
     */
    short ph;

    /**
     * Vector of {0,1,2} representing the genotypes of every variant for the
     * sample
     */
    std::vector<short> variants;

    friend std::istream &operator>>(std::istream &str, Sample &s)
    {
        std::string line;
        Sample tmp;
        if (std::getline(str, line)) {
            bool read_ok = true;
            std::stringstream iss(line);
            if (s.has_info) {
                read_ok = read_ok && iss >> tmp.fid;
                read_ok = read_ok && iss >> tmp.iid;
                read_ok = read_ok && iss >> tmp.f_iid;
                read_ok = read_ok && iss >> tmp.m_iid;
                read_ok = read_ok && iss >> tmp.sex;
            }
            for (size_t i = 0; i < s.read_mask.size(); ++i) {
                short val;
                read_ok = read_ok && iss >> val;
                if (s.read_mask[i]) {
                    if (i == s.phenotype_index) {
                        // If info is missing, phenotype is in {0,1} and needs
                        // to be added 1
                        tmp.ph = val + !s.has_info;
                    } else {
                        tmp.variants.push_back(val);
                    }
                }
            }
            if (read_ok) {
                s.swap(tmp);
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
        std::swap(variants, other.variants);
    }

  private:
    /**
     * Data came with sample information or not. If not, all other members are
     * uninitialized except for the information of the variants, and phenotypes
     * need to be adapted
     */
    bool has_info;

    /**
     * Boolean mask indicating what fields to read from the input data
     */
    std::vector<bool> read_mask;

    /**
     * Position of the phenotype value in the input information
     */
    size_t phenotype_index;
};

void read_samples(const std::string &rawfile, std::vector<Sample> &samples,
                  size_t &cases, size_t &ctrls)
{
    std::ifstream file;
    file.open(rawfile.c_str(), std::ios::in);
    if (!file.is_open()) {
        throw std::runtime_error("Error while opening " + rawfile +
                                 ", check file path/permissions");
    }
    std::string header;
    std::getline(file, header);
    // Identify fields in header
    std::regex word_regex("(\\w+)");
    auto header_it =
             std::sregex_iterator(header.begin(), header.end(), word_regex),
         header_end = std::sregex_iterator();
    // Check if first columns are sample information fields
    size_t i = 0;
    for (; header_it != header_end &&
           (header_it->str() == "FID" || header_it->str() == "IID" ||
            header_it->str() == "PAT" || header_it->str() == "MAT" ||
            header_it->str() == "SEX");
         ++i, ++header_it) {
    }
    bool has_info = i == 5;
    if (!has_info) {
        header_it =
            std::sregex_iterator(header.begin(), header.end(), word_regex);
    }
    // Create mask to select columns to read, locate phenotype column
    std::vector<bool> read_mask;
    size_t phenotype_index = -1;
    for (i = 0; header_it != header_end; ++i, ++header_it) {
        auto word = header_it->str();
        if (word == "PHENOTYPE" || word == "Class") {
            if (phenotype_index != (size_t)-1) {
                throw std::runtime_error(
                    "Error while parsing \"" + rawfile +
                    "\", multiple phenotype columns found");
            } else {
                read_mask.push_back(true);
                phenotype_index = i;
            }
        } else {
            std::smatch m;
            read_mask.push_back(std::regex_match(
                word, m, std::basic_regex("(\\w+_[ACGT])|([[:alnum:]]+)")));
        }
    }
    // Read samples
    Sample s(has_info, read_mask, phenotype_index);
    ctrls = 0;
    while (file >> s) {
        samples.push_back(s);
        ctrls += s.ph == 1;
    }
    cases = samples.size() - ctrls;
    file.close();
};

template <typename T>
void fill(std::unique_ptr<GenotypeTable<T>[]> &array,
          const std::vector<Sample> &samples, const size_t cases_words,
          const size_t ctrls_words)
{
    constexpr size_t BITS = sizeof(T) * 8; // Number of bits in T
    // Buffers
    T cases_buff[3], ctrls_buff[3];
    for (size_t i = 0; i < samples[0].variants.size(); i++) {
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
            // For each sample, check phenotype class
            if (samples[j].ph == 1) { // If it's a control append genotype to
                // the 3 control buffers
                for (auto k = 0; k < 3; k++) {
                    ctrls_buff[k] =
                        (ctrls_buff[k] << 1) + (samples[j].variants[i] == k);
                }
                ctrls_cnt++;
            } else { // Else append genotype to the 3 cases buffers
                for (auto k = 0; k < 3; k++) {
                    cases_buff[k] =
                        (cases_buff[k] << 1) + (samples[j].variants[i] == k);
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
        // Write 0 in the remaining uninitialized words of the controls
        // array
        for (auto i = (cases_cnt + BITS - 1) / BITS; i < cases_words; i++) {
            for (auto k = 0; k < 3; k++) {
                table.cases[k * cases_words + i] = 0;
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
        for (auto i = (ctrls_cnt + BITS - 1) / BITS; i < ctrls_words; i++) {
            for (auto k = 0; k < 3; k++) {
                table.ctrls[k * ctrls_words + i] = 0;
            }
        }
    }
};

template <typename T>
void read(const std::string rawfile, size_t &cases, size_t &ctrls, size_t &snps,
          std::unique_ptr<GenotypeTable<T>[]> &array)
{
    std::vector<Sample> samples;
    read_samples(rawfile, samples, cases, ctrls);
    snps = samples[0].variants.size();
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
    fill(array, samples, cases_words, ctrls_words);
}
} // namespace RawFile

#endif
