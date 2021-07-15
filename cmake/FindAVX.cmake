include(CheckCXXSourceCompiles)

set(CMAKE_REQUIRED_FLAGS_SAVE ${CMAKE_REQUIRED_FLAGS})

set(CMAKE_REQUIRED_FLAGS
    "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}")

if (NOT DEFINED AVX_ENABLED)
    set(AVX_CODE "
        #include <immintrin.h>

        int main()
        {
            __m256d a = _mm256_add_pd(_mm256_setzero_pd(),
                                      _mm256_setzero_pd());
            return 0;
        }
    ")
    check_cxx_source_compiles("${AVX_CODE}" AVX_ENABLED)
    set_property(CACHE AVX_ENABLED PROPERTY TYPE BOOL)
endif()

if (NOT DEFINED AVX2_ENABLED)
    set(AVX2_CODE "
        #include <immintrin.h>

        int main()
        {
            __m256i a = _mm256_add_epi64(_mm256_setzero_si256(),
                                         _mm256_setzero_si256());
            return 0;
        }
    ")
    check_cxx_source_compiles("${AVX2_CODE}" AVX2_ENABLED)
    set_property(CACHE AVX2_ENABLED PROPERTY TYPE BOOL)
endif()

if (NOT DEFINED FMA_ENABLED)
    set(FMA_CODE "
        #include <immintrin.h>

        int main()
        {
            __m256d a = _mm256_fmadd_pd(_mm256_setzero_pd(),
                                        _mm256_setzero_pd(),
                                        _mm256_setzero_pd());
            return 0;
        }
    ")
    check_cxx_source_compiles("${FMA_CODE}" FMA_ENABLED)
    set_property(CACHE FMA_ENABLED PROPERTY TYPE BOOL)
endif()

if (NOT DEFINED AVX512F_ENABLED)
    set(AVX512F_CODE "
        #include <immintrin.h>

        int main()
        {
            __m512d a = _mm512_add_pd(_mm512_setzero_pd(),
                                      _mm512_setzero_pd());
            return 0;
        }
    ")
    check_cxx_source_compiles("${AVX512F_CODE}" AVX512F_ENABLED)
    set_property(CACHE AVX512F_ENABLED PROPERTY TYPE BOOL)
endif()

set(CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS_SAVE})
unset(CMAKE_REQUIRED_FLAGS_SAVE)
