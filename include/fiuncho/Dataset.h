/**
 * @file Dataset.h
 * @author Christian Ponte
 * @brief Declares and implements the Dataset class
 */

#ifndef FIUNCHO_DATASET_H
#define FIUNCHO_DATASET_H

#include <filesystem>
#include <fiuncho/GenotypeTable.h>
#include <fiuncho/dataset/RAWFile.hpp>
#include <fiuncho/dataset/TPEDFile.hpp>
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

    static Dataset<T> read(const std::string file)
    {
        // Dataset inner variables
        size_t cases, ctrls, snps;
        std::unique_ptr<GenotypeTable<T>[]> array;

        const auto ext = std::filesystem::path(file).extension().string();
        if (ext == ".raw") {
            // Load data
            RawFile::read<T>(file, cases, ctrls, snps, array);
            // Create Dataset
            return Dataset<T>(cases, ctrls, snps, array);
        } else {
            throw std::runtime_error("Unrecognized input file extension " +
                                     ext);
        }
    }

    static Dataset<T> read(const std::string file1, const std::string file2)
    {
        // Declare inner variables
        size_t cases, ctrls, snps;
        std::unique_ptr<GenotypeTable<T>[]> array;

        const auto ext1 = std::filesystem::path(file1).extension().string(),
                   ext2 = std::filesystem::path(file2).extension().string();
        if (ext1 == ".tped" || ext1 == ".tfam") {
            if (ext1 == ext2 || (ext2 != ".tped" && ext2 != ".tfam")) {
                throw std::runtime_error("Expecting a tfam file after \"" +
                                         file1 + "\", but found a " + ext2);
            }
            // Load data
            if (ext1 == ".tped") {
                TPEDFile::read<T>(file1, file2, cases, ctrls, snps, array);
            } else {
                TPEDFile::read<T>(file2, file1, cases, ctrls, snps, array);
            }
            // Create Dataset
            return Dataset<T>(cases, ctrls, snps, array);
        } else {
            throw std::runtime_error("Unrecognized input file extension " +
                                     ext1);
        }
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
            std::unique_ptr<GenotypeTable<T>[]> &ptr)
        : cases(cases_count), ctrls(ctrls_count), snps(snps_count),
          data(ptr.get()), array(std::move(ptr))
    {
    }

    std::unique_ptr<GenotypeTable<T>[]> array;
};

#endif
