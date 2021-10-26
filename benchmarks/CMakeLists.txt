################################### Targets  ###################################

add_executable(bench_gtable gtable.cpp)
target_link_libraries(bench_gtable PRIVATE libfiuncho)

add_executable(bench_ctable ctable.cpp)
target_link_libraries(bench_ctable PRIVATE libfiuncho)

add_executable(bench_mi mi.cpp)
target_link_libraries(bench_mi PRIVATE libfiuncho)

add_executable(bench_threadedsearch threadedsearch.cpp)
target_link_libraries(bench_threadedsearch PRIVATE libfiuncho)
target_compile_definitions(bench_threadedsearch PRIVATE -DBENCHMARK)

add_executable(bench_nonsegmented_sequential_search
    nonsegmented_sequential_search.cpp)
target_link_libraries(bench_nonsegmented_sequential_search PRIVATE libfiuncho)

add_executable(bench_segmented_sequential_search
    segmented_sequential_search.cpp)
target_link_libraries(bench_segmented_sequential_search PRIVATE libfiuncho)
