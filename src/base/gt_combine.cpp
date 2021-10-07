#include <fiuncho/GenotypeTable.h>

template <>
void GenotypeTable<uint64_t>::combine(const GenotypeTable<uint64_t> &t1,
                                      const GenotypeTable<uint64_t> &t2,
                                      GenotypeTable<uint64_t> &out) noexcept
{
    size_t i, j, k;
    // Constant references to objects are not enough for the compiler to
    // determine the number of loop iterations for vectorization, therefore they
    // have to be copied to a local constant
    const size_t table_size = t1.size;
    const size_t cases_words = t1.cases_words;
    const size_t ctrls_words = t1.ctrls_words;
    // Compute genotype tables for cases
    for (i = 0; i < table_size; i++) {
        for (j = 0; j < 3; j++) {
            for (k = 0; k < cases_words; k++) {
                out.cases[(i * 3 + j) * cases_words + k] =
                    t1.cases[i * cases_words + k] &
                    t2.cases[j * cases_words + k];
            }
        }
    }
    // Compute genotype tables for ctrls
    for (i = 0; i < table_size; i++) {
        for (j = 0; j < 3; j++) {
            for (k = 0; k < ctrls_words; k++) {
                out.ctrls[(i * 3 + j) * ctrls_words + k] =
                    t1.ctrls[i * ctrls_words + k] &
                    t2.ctrls[j * ctrls_words + k];
            }
        }
    }
}
