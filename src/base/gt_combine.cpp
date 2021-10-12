#include <fiuncho/GenotypeTable.h>

inline void combine_aux(const uint64_t *gt_tbl1, const size_t gt_size,
                        const uint64_t *gt_tbl2, const size_t words,
                        uint64_t *gt_tbl3)
{
    size_t i, j, k;
    // Combine two genotype tables and store it in a larger table
    // Compute genotype tables for cases
    for (i = 0; i < gt_size; i++) {
        for (j = 0; j < 3; j++) {
            for (k = 0; k < words; k++) {
                gt_tbl3[(i * 3 + j) * words + k] =
                    gt_tbl1[i * words + k] &
                    gt_tbl2[j * words + k];
            }
        }
    }
}

template <>
void GenotypeTable<uint64_t>::combine(const GenotypeTable<uint64_t> &t1,
                                      const GenotypeTable<uint64_t> &t2,
                                      GenotypeTable<uint64_t> &out) noexcept
{
    combine_aux(t1.cases, t1.size, t2.cases, t1.cases_words, out.cases);
    combine_aux(t1.ctrls, t1.size, t2.ctrls, t1.ctrls_words, out.ctrls);
}
