#include <gtest/gtest.h>
#include <fiuncho/ContingencyTable.h>
#include <fiuncho/algorithms/MutualInformation.h>

namespace
{
TEST(MI, compute)
{
    ContingencyTable<uint32_t> ctable(2);
#ifdef ALIGN
    for (size_t i = 9; i < ctable.size; i++) {
        ctable.cases[i] = 0;
        ctable.ctrls[i] = 0;
    }
#endif
    MutualInformation<float> mi(9, 18);

    for (size_t i = 0; i < 9; i++) {
        ctable.cases[i] = 1;
        ctable.ctrls[i] = 2;
    }
    EXPECT_NEAR(0, mi.compute(ctable), 1E-5);

    for (size_t i = 0; i < 9; i++) {
        ctable.cases[i] = 0;
        ctable.ctrls[i] = 0;
    }
    ctable.cases[0] = 9;
    ctable.ctrls[1] = 18;
    EXPECT_NEAR(0.6365141682948129, mi.compute(ctable), 1E-5);
}
} // namespace
