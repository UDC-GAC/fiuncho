#include <iostream>
#include <string>
#include <thread>
#include <time.h>
#include <vector>

/**
 * @brief Meassure elapsed time of a generic function pinned to a particular CPU
 * core.
 *
 * @tparam T Function type (should be deduced by the compiler)
 * @tparam Args Arguments to the function T (should also be deduced)
 * @param affinity CPU core to which the thread should be pinned
 * @param barrier Pthread barrier used to synchronize the different threads
 * @param[out] time Time elapsed during function func
 * @param func Function to meassure
 * @param args Arguments to the function func
 */

template <typename T, typename... Args>
void pinned_time(const int affinity, pthread_barrier_t &barrier, double &time,
                 T func, Args &&...args)
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
    // Write time into output variable
    time =
        end.tv_sec + end.tv_nsec * 1E-9 - start.tv_sec - start.tv_nsec * 1E-9;
}

/**
 * @brief Run a multithreaded benchmark, measuring the elapsed time of each
 * thread and pinning the threads to a particular CPU affinity.
 *
 * @tparam @tparam T Function type (should be deduced by the compiler)
 * @tparam Args Arguments to the function T (should also be deduced)
 * @param affinity Vector of ints containing the CPU affinity to be used
 * @param func Function to meassure
 * @param args Arguments to the function func
 * @return std::vector<double> Elapsed time of the different threads spent in
 * function func
 */

template <typename T, typename... Args>
std::vector<double> multithread_pinned_time(const std::vector<int> &affinity,
                                            T func, Args &&...args)
{
    const int thread_count = affinity.size();
    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier, NULL, thread_count);
    std::vector<std::thread> threads;
    std::vector<double> times(thread_count);
    // Launch thread_count-1 threads
    for (auto i = 1; i < thread_count; ++i) {
        threads.emplace_back([&, i]() {
            pinned_time(affinity[i], barrier, times[i], func,
                        std::forward<Args>(args)...);
        });
    }
    // Use current thread as well
    pinned_time(affinity[0], barrier, times[0], func,
                std::forward<Args>(args)...);
    // Wait for all threads to complete
    for (auto &t : threads) {
        t.join();
    }
    return times;
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
