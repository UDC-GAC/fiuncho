/*
 * MPI3SNP ~ https://github.com/chponte/mpi3snp
 *
 * Copyright 2018 Christian Ponte
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
 * documentation files (the "Software"), to deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the
 * Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
 * WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * @file Statistics.h
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief Statistics class definition and private template functions implementation.
 */

#ifndef MPI3SNP_STATISTICS_H
#define MPI3SNP_STATISTICS_H

#include <vector>
#include <tuple>
#include <string>
#include <pthread.h>

class Statistics {
public:
    Statistics();

    ~Statistics();

    inline void Addi(const std::string &label, int value) {
        pthread_mutex_lock(&ints_mutex);
        Add<int>(label, value, ints);
        pthread_mutex_unlock(&ints_mutex);
    }

    inline int Geti(const std::string &label) {
        int retval = 0;
        pthread_mutex_lock(&ints_mutex);
        retval = Get<int>(label, ints);
        pthread_mutex_unlock(&ints_mutex);
        return retval;
    }

    inline void Addl(const std::string &label, long value) {
        pthread_mutex_lock(&longs_mutex);
        Add<long>(label, value, longs);
        pthread_mutex_unlock(&longs_mutex);
    }

    inline long Getl(const std::string &label) {
        long retval = 0;
        pthread_mutex_lock(&longs_mutex);
        retval = Get<long>(label, longs);
        pthread_mutex_unlock(&longs_mutex);
        return retval;
    }

    void Begin_timer(const std::string &label);

    double End_timer(const std::string &label);

    double Get_timer(const std::string &label);

    std::string To_string();

private:
    template<typename T>
    static typename T::const_iterator Find(const std::string &label, const T &vector) {
        auto it = vector.begin();
        while (it < vector.end() && std::get<0>(*it).compare(label) != 0) {
            it++;
        }
        return it;
    }

    template<typename T>
    static void Add(const std::string &label, const T &value, std::vector<std::pair<std::string, T>> &vector) {
        auto pos = Find<std::vector<std::pair<std::string, T>>>(label, vector);
        if (pos != vector.end()) {
            pos = vector.erase(pos);
        }
        vector.insert(pos, std::make_pair(label, value));
    }

    template<typename T>
    static T Get(const std::string &label, const std::vector<std::pair<std::string, T>> &vector) {
        auto pos = Find<std::vector<std::pair<std::string, T>>>(label, vector);
        if (pos != vector.end()) {
            return std::get<1>(*pos);
        }
        return 0;
    }

    pthread_mutex_t ints_mutex;
    std::vector<std::pair<std::string, int>> ints;

    pthread_mutex_t longs_mutex;
    std::vector<std::pair<std::string, long>> longs;

    std::vector<std::tuple<std::string, double, bool>>::const_iterator Find_timer_label(const std::string &label);

    pthread_mutex_t timers_mutex;
    std::vector<std::tuple<std::string, double, bool>> timers;
};

#endif //MPI3SNP_STATISTICS_H
