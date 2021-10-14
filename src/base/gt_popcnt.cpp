#include <bitset>
#include <fiuncho/GenotypeTable.h>

inline void popcnt_aux(const uint64_t *gt_tbl1, const size_t gt_size,
                       const uint64_t *gt_tbl2, const size_t words,
                       uint32_t *ct_tbl, const size_t ct_size)
{
    size_t i, j, k;
    // Compute contingency tables for cases
    for (i = 0; i < gt_size; i++) {
        for (j = 0; j < 3; j++) {
            ct_tbl[i * 3 + j] = 0;
            for (k = 0; k < words; k++) {
                ct_tbl[i * 3 + j] += std::bitset<64>(gt_tbl1[i * words + k] &
                                                     gt_tbl2[j * words + k])
                                         .count();
            }
        }
    }
    for (i = gt_size * 3; i < ct_size; ++i) {
        ct_tbl[i] = 0;
    }
}

template <>
template <>
void GenotypeTable<uint64_t>::combine_and_popcnt(
    const GenotypeTable<uint64_t> &t1, const GenotypeTable<uint64_t> &t2,
    ContingencyTable<uint32_t> &out) noexcept
{
    popcnt_aux(t1.cases, t1.size, t2.cases, t1.cases_words, out.cases,
               out.size);
    popcnt_aux(t1.ctrls, t1.size, t2.ctrls, t1.ctrls_words, out.ctrls,
               out.size);
}
