#ifndef FIUNCHO_ALGORITHM_H
#define FIUNCHO_ALGORITHM_H

#include <fiuncho/ContingencyTable.h>

template <class T>
class Algorithm {
  public:
    virtual ~Algorithm() = default;

    template <class U>
    T compute(const ContingencyTable<U> &table) const noexcept;
};

#endif
