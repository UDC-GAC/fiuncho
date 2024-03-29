/**
 * @file Dataset.h
 * @date 01/11/2021
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
     * @param inputs List of files to read (in any order). Supported files are:
     *      1. <code>tped</code> + <code>tfam</code> files
     *      2. <code>raw</code> files
     * @return A Dataset object
     */

    static Dataset<T> read(const std::vector<std::string> &inputs)
    {
        // Declare inner variables
        size_t cases, ctrls, snps;
        std::unique_ptr<GenotypeTable<T>[]> array;

        switch (inputs.size()) {
        case 1: {
            const auto file = inputs[0];
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
        case 2: {
            const auto file1 = inputs[0], file2 = inputs[1];
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
        default:
            throw std::runtime_error("Unrecognized number of input files to "
                                     "Dataset factory constructor");
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
