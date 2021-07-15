include(CheckCXXSourceCompiles)

set(CMAKE_REQUIRED_FLAGS_SAVE ${CMAKE_REQUIRED_FLAGS})

set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}")

if (NOT DEFINED SVML_AVAILABLE)
    set(SVML_CODE "
        #include <immintrin.h>

        int main()
        {
            __m256d a = _mm256_log_pd(_mm256_setzero_pd());
            return 0;
        }
    ")
    check_cxx_source_compiles("${SVML_CODE}" SVML_AVAILABLE)
    set_property(CACHE SVML_AVAILABLE PROPERTY TYPE BOOL)
endif()

if (NOT DEFINED LIBMVEC_AVAILABLE)
    set(LIBMVEC_CODE "
        #include <immintrin.h>

        extern \"C\" {
            __m256d _ZGVdN4v_log(__m256d x);
        }

        int main()
        {
            __m256d a = _ZGVdN4v_log(_mm256_setzero_pd());
            return 0;
        }
    ")
    check_cxx_source_compiles("${LIBMVEC_CODE}" LIBMVEC_AVAILABLE)
    set_property(CACHE LIBMVEC_AVAILABLE PROPERTY TYPE BOOL)
endif()

set(CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS_SAVE})
unset(CMAKE_REQUIRED_FLAGS_SAVE)
