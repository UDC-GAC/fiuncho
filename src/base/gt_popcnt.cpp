#include <bitset>
#include <fiuncho/GenotypeTable.h>

template <>
template <>
void GenotypeTable<uint64_t>::combine_and_popcnt(
    const GenotypeTable<uint64_t> &t1, const GenotypeTable<uint64_t> &t2,
    ContingencyTable<uint32_t> &out) noexcept
{
    size_t i, j, k;
    // Constant references to objects are not enough for the compiler to
    // determine the number of loop iterations for vectorization, therefore they
    // have to be copied to a local constant
    const size_t table_size = t1.size;
    const size_t out_size = out.size;
    const size_t cases_words = t1.cases_words;
    const size_t ctrls_words = t1.ctrls_words;
    // Set all contingency table values to 0
    for (i = 0; i < out_size; i++) {
        out.cases[i] = 0;
        out.ctrls[i] = 0;
    }
    // Compute contingency tables for cases
    for (i = 0; i < table_size; i++) {
        for (j = 0; j < 3; j++) {
            for (k = 0; k < cases_words; k++) {
                out.cases[i * 3 + j] +=
                    std::bitset<64>(t1.cases[i * cases_words + k] &
                                    t2.cases[j * cases_words + k])
                        .count();
            }
        }
    }
    // Compute contingency tables for ctrls
    for (i = 0; i < table_size; i++) {
        for (j = 0; j < 3; j++) {
            for (k = 0; k < ctrls_words; k++) {
                out.ctrls[i * 3 + j] +=
                    std::bitset<64>(t1.ctrls[i * ctrls_words + k] &
                                    t2.ctrls[j * ctrls_words + k])
                        .count();
            }
        }
    }
}
