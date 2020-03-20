/*
 * This file is part of Fiuncho.
 *
 * Fiuncho is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Fiuncho is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Fiuncho. If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * @file Statistics.h
 * @author Christian Ponte
 * @date 28 May 2018
 *
 * @brief Statistics class definition and private template functions
 * implementation.
 */

#ifndef FIUNCHO_STATISTICS_H
#define FIUNCHO_STATISTICS_H

#include <pthread.h>
#include <string>
#include <tuple>
#include <vector>

class Statistics {
  public:
    Statistics();

    ~Statistics();

    inline void addi(const std::string &label, int value) {
        pthread_mutex_lock(&ints_mutex);
        add<int>(label, value, ints);
        pthread_mutex_unlock(&ints_mutex);
    }

    inline int geti(const std::string &label) {
        int retval = 0;
        pthread_mutex_lock(&ints_mutex);
        retval = get<int>(label, ints);
        pthread_mutex_unlock(&ints_mutex);
        return retval;
    }

    inline void addl(const std::string &label, long value) {
        pthread_mutex_lock(&longs_mutex);
        add<long>(label, value, longs);
        pthread_mutex_unlock(&longs_mutex);
    }

    inline long getl(const std::string &label) {
        long retval = 0;
        pthread_mutex_lock(&longs_mutex);
        retval = get<long>(label, longs);
        pthread_mutex_unlock(&longs_mutex);
        return retval;
    }

    void begin_timer(const std::string &label);

    double end_timer(const std::string &label);

    double get_timer(const std::string &label);

    std::string str();

  private:
    template <typename T>
    static typename T::const_iterator find(const std::string &label,
                                           const T &vector) {
        auto it = vector.begin();
        while (it < vector.end() && std::get<0>(*it) != label) {
            it++;
        }
        return it;
    }

    template <typename T>
    static void add(const std::string &label, const T &value,
                    std::vector<std::pair<std::string, T>> &vector) {
        auto pos = find<std::vector<std::pair<std::string, T>>>(label, vector);
        if (pos != vector.end()) {
            pos = vector.erase(pos);
        }
        vector.insert(pos, std::make_pair(label, value));
    }

    template <typename T>
    static T get(const std::string &label,
                 const std::vector<std::pair<std::string, T>> &vector) {
        auto pos = find<std::vector<std::pair<std::string, T>>>(label, vector);
        if (pos != vector.end()) {
            return std::get<1>(*pos);
        }
        return 0;
    }

    pthread_mutex_t ints_mutex;
    std::vector<std::pair<std::string, int>> ints;

    pthread_mutex_t longs_mutex;
    std::vector<std::pair<std::string, long>> longs;

    std::vector<std::tuple<std::string, double, bool>>::const_iterator
    find_timer_label(const std::string &label);

    pthread_mutex_t timers_mutex;
    std::vector<std::tuple<std::string, double, bool>> timers;
};

#endif
