include(CheckCXXSourceCompiles)

set(CMAKE_REQUIRED_FLAGS_SAVE ${CMAKE_REQUIRED_FLAGS})
set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS} -g")

set(TEMPLATE_CODE "
#include <type_traits>

int main()
{
    static_assert(
#ifdef %DEFINITION%
        true,
#else
        false,
#endif
        \"flag not enabled\");
    return 0;
}")
set(CPU_EXTENSIONS
    "FMA"
    "FMA4"
    "AVX"
    "AVX2"
    "AVX512F"
    "AVX512DQ"
    "AVX512IFMA"
    "AVX512CD"
    "AVX512BW"
    "AVX512VL"
    "AVX512VBMI"
    "AVX512VBMI2"
    "AVX512VNNI"
    "AVX512VNNIW"
    "AVX512BITALG"
    "AVX5124FMAPS"
    "AVX512VPOPCNTDQ"
    "AVX512VP2INTERSECT")

foreach(E IN LISTS CPU_EXTENSIONS)
    string(REPLACE "%DEFINITION%" "__${E}__" ${E}_CODE "${TEMPLATE_CODE}")
    check_cxx_source_compiles("${${E}_CODE}" ${E}_ENABLED)
    set_property(CACHE ${E}_ENABLED PROPERTY TYPE BOOL)
endforeach()

set(CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS_SAVE})
unset(CMAKE_REQUIRED_FLAGS_SAVE)
