#include <iostream>
#include <string>
#include <thread>
#include <time.h>
#include <vector>

pthread_barrier_t barrier;

/**
 * @brief Meassure elapsed time of a generic function pinned to a particular CPU
 * core.
 *
 * @tparam T Function type (should be deduced by the compiler)
 * @tparam Args Arguments to the function T (should also be deduced)
 * @param affinity CPU core to which the thread should be pinned
 * @param func Function to meassure
 * @param args Arguments to the function func
 * @return double Time elapsed during function func
 */

template <typename T, typename... Args>
double pinned_time(const int affinity, T func, Args &&...args)
{
    // Pin thread to CPU core
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(affinity, &cpuset);
    int rc = pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
    if (rc != 0) {
        std::cerr << "Error calling pthread_setaffinity_np: " << rc << "\n";
    }
    // Sync threads
    pthread_barrier_wait(&barrier);
    // Measure search time
    struct timespec start, end;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start);
    func(std::forward<Args>(args)...);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end);
    return end.tv_sec + end.tv_nsec * 1E-9 - start.tv_sec -
           start.tv_nsec * 1E-9;
}

std::vector<int> split_into_ints(const std::string &s, const char sep)
{
    std::vector<int> ints;
    size_t pos = s.find(sep, 0), prev = 0;
    while (pos != std::string::npos) {
        ints.push_back(atoi(s.substr(prev, pos - prev).c_str()));
        prev = pos + 1;
        pos = s.find(sep, pos + 1);
    }
    ints.push_back(atoi(s.substr(prev).c_str()));
    return ints;
}
