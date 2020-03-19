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
 * @file Statistics.cpp
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief Statistics class members implementation.
 */

#include "Statistics.h"
#include <mpi.h>

Statistics::Statistics() {
    pthread_mutex_init(&ints_mutex, nullptr);
    pthread_mutex_init(&longs_mutex, nullptr);
    pthread_mutex_init(&timers_mutex, nullptr);
}

Statistics::~Statistics() {
    pthread_mutex_destroy(&ints_mutex);
    pthread_mutex_destroy(&longs_mutex);
    pthread_mutex_destroy(&timers_mutex);
}

void Statistics::Begin_timer(const std::string &label) {
    pthread_mutex_lock(&timers_mutex);
    auto pos = Find_timer_label(label);
    if (pos != timers.end()){
        pos = timers.erase(pos);
    }
    timers.insert(pos, std::make_tuple(label, MPI_Wtime(), true));
    pthread_mutex_unlock(&timers_mutex);
}

// TODO: raise exception if label does not exist
double Statistics::End_timer(const std::string &label) {
    pthread_mutex_lock(&timers_mutex);
    auto pos = Find_timer_label(label);
    if (pos != timers.end()) {
        double time;
        bool active;
        std::tie(std::ignore, time, active) = *pos;

        if (active) {
            time = MPI_Wtime() - time;
            pos = timers.erase(pos);
            timers.insert(pos, std::make_tuple(label, time, false));
            pthread_mutex_unlock(&timers_mutex);
            return time;
        }
    }
    pthread_mutex_unlock(&timers_mutex);
    return 0;
}

// TODO: raise exception if label does not exist
double Statistics::Get_timer(const std::string &label) {
    double out = 0;
    pthread_mutex_lock(&timers_mutex);
    auto pos = Find_timer_label(label);
    if (pos != timers.end()) {
        out = std::get<1>(*pos);
    }
    pthread_mutex_unlock(&timers_mutex);
    return out;
}

std::string Statistics::To_string() {
    std::string output("Statistics\n");
    for (auto item : timers){
        if (!std::get<2>(item)) {
            output += "\t" + std::get<0>(item) + ": " + std::to_string(std::get<1>(item)) + " seconds\n";
        }
    }

    for (auto item : ints){
        output += "\t" + item.first + ": " + std::to_string(item.second) + "\n";
    }
    for (auto item : longs){
        output += "\t" + item.first + ": " + std::to_string(item.second) + "\n";
    }

    return output;
}

std::vector<std::tuple<std::string, double, bool>>::const_iterator
Statistics::Find_timer_label(const std::string &label) {
    auto it = timers.begin();
    while (it < timers.end() && std::get<0>(*it).compare(label) != 0) {
        it++;
    }
    return it;
}