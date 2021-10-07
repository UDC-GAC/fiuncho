#include <cmath>
#include <fiuncho/GenotypeTable.h>

template <>
GenotypeTable<uint64_t>::GenotypeTable(const short order,
                                       const size_t cases_words,
                                       const size_t ctrls_words)
    : order(order), size(std::pow(3, order)), cases_words(cases_words),
      ctrls_words(ctrls_words),
      alloc(std::make_unique<uint64_t[]>(size * (cases_words + ctrls_words))),
      cases(alloc.get()), ctrls(cases + size * cases_words)
{
}
