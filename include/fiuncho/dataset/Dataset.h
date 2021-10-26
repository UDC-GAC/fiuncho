/**
 * @file Dataset.h
 * @author Christian Ponte
 * @brief Declares and implements the Dataset class
 */

#ifndef FIUNCHO_DATASET_H
#define FIUNCHO_DATASET_H

#include <array>
#include <fiuncho/GenotypeTable.h>
#include <fiuncho/dataset/Individual.h>
#include <fiuncho/dataset/SNP.h>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

/**
 * @class Dataset
 * @brief Class representing a collection of \a N GenotypeTable's, each
 * containing a single SNP, for \a M individuals
 *
 * @tparam T data type used to represent the individual information in the
 * GenotypeTable's
 */

template <class T> class Dataset
{
  public:
    Dataset(const Dataset<T> &) = delete;
    Dataset(Dataset<T> &&) = default;

    /**
     * @name Factory methods
     */
    //@{

    /**
     * Read input data and store it using a GenotypeTable representation. The
     * underlying arrays used in the different tables are allocated contiguously
     * in memory.
     *
     * @param tped Path to the tped input file
     * @param tfam Path to the tfam input file
     * @return A Dataset object
     */

    static Dataset<T> read(std::string tped, std::string tfam)
    {
        std::vector<Individual> individuals;
        std::vector<SNP> snps;
        size_t cases_count, ctrls_count;
        // Read input files
        read_individuals(tfam, individuals, cases_count, ctrls_count);
        read_snps(tped, individuals, snps);
        // Allocate enough space for representing all SNPs for all individuals
#ifdef ALIGN
        constexpr size_t NT = ALIGN / sizeof(T); // Number of T's in ALIGN bytes
        constexpr size_t NBITS = ALIGN * 8; // Number of bits in ALIGN bytes
        const size_t cases_words = (cases_count + NBITS - 1) / NBITS * NT,
                     ctrls_words = (ctrls_count + NBITS - 1) / NBITS * NT;
#else
        // Allocate enough space for representing all SNPs for all individuals
        constexpr size_t NBITS = sizeof(T) * 8; // Number of bits in T
        const size_t cases_words = (cases_count + NBITS - 1) / NBITS,
                     ctrls_words = (ctrls_count + NBITS - 1) / NBITS;
#endif
        // Create Dataset object
        Dataset<T> d(cases_count, ctrls_count, snps.size(),
                     GenotypeTable<T>::make_array(snps.size(), 1, cases_words,
                                                  ctrls_words));
        // Fill it with the information of individuals and snps
        populate(d, individuals, snps, cases_words, ctrls_words);
        return d;
    }

    //@}

    /**
     * @name Methods
     */
    //@{

    /**
     * Access the bit table vector.
     *
     * @param i Index of the table in the vector
     * @return A reference to the GenotypeTable object
     */
    GenotypeTable<T> &operator[](int i) { return data[i]; }

    const GenotypeTable<T> &operator[](int i) const { return data[i]; }

    //@}

    /**
     * @name Attributes
     */
    //@{

    /**
     * Number of individuals in the case group of the data set
     */
    const size_t cases;

    /**
     * Number of individuals in the control group of the data set
     */
    const size_t ctrls;

    /**
     * Number of SNPs in the data set
     */
    const size_t snps;

    /**
     * Vector holding all the GenotypeTable's representing the individual SNP's
     * information
     */
    GenotypeTable<T> *data;

    //@}

  private:
    Dataset(size_t cases_count, size_t ctrls_count, size_t snps_count,
            std::unique_ptr<GenotypeTable<T>[]> ptr)
        : cases(cases_count), ctrls(ctrls_count), snps(snps_count),
          data(ptr.get()), array(std::move(ptr))
    {
    }

    inline static void read_individuals(const std::string &tfam,
                                        std::vector<Individual> &individuals,
                                        size_t &cases, size_t &ctrls)
    {
        std::ifstream file;
        file.open(tfam.c_str(), std::ios::in);
        if (!file.is_open()) {
            throw std::runtime_error("Error while opening " + tfam +
                                     ", check file path/permissions");
        }
        try {
            Individual ind;
            ctrls = 0;
            while (file >> ind) {
                individuals.push_back(ind);
                ctrls += ind.ph == 1;
            }
            cases = individuals.size() - ctrls;
        } catch (const Individual::InvalidIndividual &e) {
            throw std::runtime_error("Error in " + tfam + ":" +
                                     std::to_string(individuals.size() + 1) +
                                     ": " + e.what());
        }
        file.close();
    }

    inline static void read_snps(const std::string &tped,
                                 const std::vector<Individual> &individuals,
                                 std::vector<SNP> &snps)
    {
        std::ifstream file;
        file.open(tped.c_str(), std::ios::in);
        if (!file.is_open()) {
            throw std::runtime_error("Error while opening " + tped +
                                     ", check file path/permissions");
        }
        try {
            SNP snp;
            while (file >> snp) {
                if (snp.genotypes.size() == individuals.size()) {
                    snps.push_back(snp);
                } else {
                    throw std::runtime_error(
                        "Error in " + tped + ":" +
                        std::to_string(snps.size() + 1) +
                        ": the number of nucleotides does not match "
                        "the number of individuals");
                }
            }
        } catch (const SNP::InvalidSNP &e) {
            throw std::runtime_error("Error in " + tped + ":" +
                                     std::to_string(snps.size() + 1) + ": " +
                                     e.what());
        }
        file.close();
    }

    inline static void populate(Dataset<T> &dataset,
                                const std::vector<Individual> &inds,
                                const std::vector<SNP> &snps,
                                const size_t cases_words,
                                const size_t ctrls_words)
    {
        constexpr size_t BITS = sizeof(T) * 8; // Number of bits in T
        // Buffers
        T cases_buff[3], ctrls_buff[3];
        for (size_t i = 0; i < dataset.snps; i++) {
            // Clear buffers
            for (auto k = 0; k < 3; k++) {
                cases_buff[k] = 0;
                ctrls_buff[k] = 0;
            }
            // Obtain reference to the next genotype table
            auto &table = dataset[i];
            // Populate bit table with the snp information
            size_t cases_cnt = 0;
            size_t ctrls_cnt = 0;
            for (size_t j = 0; j < inds.size(); j++) {
                // For each individual, check phenotype class
                if (inds[j].ph == 1) { // If it's a control append genotype to
                    // the 3 control buffers
                    for (auto k = 0; k < 3; k++) {
                        ctrls_buff[k] =
                            (ctrls_buff[k] << 1) + (snps[i].genotypes[j] == k);
                    }
                    ctrls_cnt++;
                } else { // Else append genotype to the 3 cases buffers
                    for (auto k = 0; k < 3; k++) {
                        cases_buff[k] =
                            (cases_buff[k] << 1) + (snps[i].genotypes[j] == k);
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
    }

    std::unique_ptr<GenotypeTable<T>[]> array;
};

#endif
