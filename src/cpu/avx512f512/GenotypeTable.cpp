#include <cmath>
#include <fiuncho/GenotypeTable.h>
#include <immintrin.h>

template <>
GenotypeTable<uint64_t>::GenotypeTable(const short order,
                                       const size_t cases_words,
                                       const size_t ctrls_words)
    : order(order), size(std::pow(3, order)), cases_words(cases_words),
      ctrls_words(ctrls_words), alloc(std::make_unique<uint64_t[]>(
                                    size * (cases_words + ctrls_words) + 8)),
      cases((uint64_t *)((((uintptr_t)alloc.get()) + 63) / 64 * 64)),
      ctrls(cases + size * cases_words)
{
}

template <>
void GenotypeTable<uint64_t>::combine(const GenotypeTable<uint64_t> &t1,
                                      const GenotypeTable<uint64_t> &t2,
                                      GenotypeTable<uint64_t> &out) noexcept
{
    size_t i, j, k;
    // Compute bit tables for cases
    for (i = 0; i < t1.size; i += 3) {
        for (j = 0; j < 3; j++) {
            for (k = 0; k < t1.cases_words; k += 8) {
                __m512i z0 = _mm512_load_si512(
                    (__m512i *)(t2.cases + j * t1.cases_words + k));
                __m512i z1 = _mm512_load_si512(
                    (__m512i *)(t1.cases + (i + 0) * t1.cases_words + k));
                __m512i z2 = _mm512_load_si512(
                    (__m512i *)(t1.cases + (i + 1) * t1.cases_words + k));
                __m512i z3 = _mm512_load_si512(
                    (__m512i *)(t1.cases + (i + 2) * t1.cases_words + k));
                __m512i z4 = _mm512_and_si512(z0, z1);
                __m512i z5 = _mm512_and_si512(z0, z2);
                __m512i z6 = _mm512_and_si512(z0, z3);
                _mm512_store_si512(
                    (__m512i *)(out.cases + ((i + j) * 3 + 0) * t1.cases_words +
                                k),
                    z4);
                _mm512_store_si512(
                    (__m512i *)(out.cases + ((i + j) * 3 + 1) * t1.cases_words +
                                k),
                    z5);
                _mm512_store_si512(
                    (__m512i *)(out.cases + ((i + j) * 3 + 2) * t1.cases_words +
                                k),
                    z6);
            }
        }
    }
    // Compute bit tables for ctrls
    for (i = 0; i < t1.size; i += 3) {
        for (j = 0; j < 3; j++) {
            for (k = 0; k < t1.ctrls_words; k += 8) {
                __m512i z0 = _mm512_load_si512(
                    (__m512i *)(t2.ctrls + j * t1.ctrls_words + k));
                __m512i z1 = _mm512_load_si512(
                    (__m512i *)(t1.ctrls + (i + 0) * t1.ctrls_words + k));
                __m512i z2 = _mm512_load_si512(
                    (__m512i *)(t1.ctrls + (i + 1) * t1.ctrls_words + k));
                __m512i z3 = _mm512_load_si512(
                    (__m512i *)(t1.ctrls + (i + 2) * t1.ctrls_words + k));
                __m512i z4 = _mm512_and_si512(z0, z1);
                __m512i z5 = _mm512_and_si512(z0, z2);
                __m512i z6 = _mm512_and_si512(z0, z3);
                _mm512_store_si512(
                    (__m512i *)(out.ctrls + ((i + j) * 3 + 0) * t1.ctrls_words +
                                k),
                    z4);
                _mm512_store_si512(
                    (__m512i *)(out.ctrls + ((i + j) * 3 + 1) * t1.ctrls_words +
                                k),
                    z5);
                _mm512_store_si512(
                    (__m512i *)(out.ctrls + ((i + j) * 3 + 2) * t1.ctrls_words +
                                k),
                    z6);
            }
        }
    }
}

template <>
template <>
void GenotypeTable<uint64_t>::combine_and_popcnt(
    const GenotypeTable<uint64_t> &t1, const GenotypeTable<uint64_t> &t2,
    ContingencyTable<uint32_t> &out) noexcept
{
    const __m512i lookup = _mm512_setr_epi64(
        0x0302020102010100llu, 0x0403030203020201llu, 0x0302020102010100llu,
        0x0403030203020201llu, 0x0302020102010100llu, 0x0403030203020201llu,
        0x0302020102010100llu, 0x0403030203020201llu);
    const __m512i low_mask = _mm512_set1_epi8(0x0f);
    size_t i, j, k, l;
    // Compute count tables for cases
    for (i = 0; i < t1.size; ++i) {
        for (j = 0; j < 3; ++j) {
            k = 0;
            __m512i acc = _mm512_setzero_si512();
            while (k < t1.cases_words) {
                __m512i local = _mm512_setzero_si512();
                for (l = 0; l < 255 / 8 && k < t1.cases_words; ++l, k += 8) {
                    __m512i z0 = _mm512_load_si512(
                        (__m512i *)(t2.cases + j * t1.cases_words + k));
                    __m512i z1 = _mm512_load_si512(
                        (__m512i *)(t1.cases + i * t1.cases_words + k));
                    __m512i z4 = _mm512_and_si512(z0, z1);
                    const __m512i lo = _mm512_and_si512(z4, low_mask);
                    const __m512i hi =
                        _mm512_and_si512(_mm512_srli_epi32(z4, 4), low_mask);
                    const __m512i popcnt1 = _mm512_shuffle_epi8(lookup, lo);
                    const __m512i popcnt2 = _mm512_shuffle_epi8(lookup, hi);
                    local = _mm512_add_epi8(local, popcnt1);
                    local = _mm512_add_epi8(local, popcnt2);
                }
                acc = _mm512_add_epi64(
                    acc, _mm512_sad_epu8(local, _mm512_setzero_si512()));
            }
            out.cases[i * 3 + j] = _mm512_reduce_add_epi64(acc);
        }
    }
    for (i = t1.size * 3; i < out.size; ++i) {
        out.cases[i] = 0;
    }
    // Compute count tables for ctrls
    for (i = 0; i < t1.size; ++i) {
        for (j = 0; j < 3; ++j) {
            k = 0;
            __m512i acc = _mm512_setzero_si512();
            while (k < t1.ctrls_words) {
                __m512i local = _mm512_setzero_si512();
                for (l = 0; l < 255 / 8 && k < t1.ctrls_words; ++l, k += 8) {
                    __m512i z0 = _mm512_load_si512(
                        (__m512i *)(t2.ctrls + j * t1.ctrls_words + k));
                    __m512i z1 = _mm512_load_si512(
                        (__m512i *)(t1.ctrls + i * t1.ctrls_words + k));
                    __m512i z4 = _mm512_and_si512(z0, z1);
                    const __m512i lo = _mm512_and_si512(z4, low_mask);
                    const __m512i hi =
                        _mm512_and_si512(_mm512_srli_epi32(z4, 4), low_mask);
                    const __m512i popcnt1 = _mm512_shuffle_epi8(lookup, lo);
                    const __m512i popcnt2 = _mm512_shuffle_epi8(lookup, hi);
                    local = _mm512_add_epi8(local, popcnt1);
                    local = _mm512_add_epi8(local, popcnt2);
                }
                acc = _mm512_add_epi64(
                    acc, _mm512_sad_epu8(local, _mm512_setzero_si512()));
            }
            out.ctrls[i * 3 + j] = _mm512_reduce_add_epi64(acc);
        }
    }
    for (i = t1.size * 3; i < out.size; ++i) {
        out.ctrls[i] = 0;
    }
}
