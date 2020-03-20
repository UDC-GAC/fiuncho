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
 * @file IOMpi.h
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief IOMpi class definition.
 */

#ifndef FIUNCHO_IOMPI_H
#define FIUNCHO_IOMPI_H

#include <cstdarg>
#include <iostream>
#include <pthread.h>

/**
 * @class IOMpi
 *
 * @brief Handles all the I/O throughout the program. Uses Meyers' Singleton
 * pattern, which describes a thread-safe singleton initialization (C++11 or
 * superior). The object destruction is handled at the end of the program.
 */

class IOMpi {
  public:
    enum Level {
        E = 0, // Error
        N = 1, // None
        B = 2, // Benchmarking
        D = 3  // Debug
    };

    static IOMpi &Instance() {
        static IOMpi inst;
        return inst;
    }

    static void Set_print_level(Level l) { Instance().level = l; }

    // delete copy and move constructors and assign operators
    IOMpi(IOMpi const &) = delete;            // Copy construct
    IOMpi(IOMpi &&) = delete;                 // Move construct
    IOMpi &operator=(IOMpi const &) = delete; // Copy assign
    IOMpi &operator=(IOMpi &&) = delete;      // Move assign

    template <Level l> inline void print(const std::string &s) {
        if (l <= level) {
            sprint_nolr(std::cout, s);
        }
    }

    template <Level l> inline void mprint(const std::string &s) {
        if (l <= level && my_rank == io_rank) {
            sprint_nolr(std::cout, s);
        }
    }

    template <Level l> inline void cprint(const std::string &s) {
        if (l <= level) {
            scprint_nol(std::cout, s);
        }
    }

    template <Level l>
    inline void sprint(std::ostream &ostream, const std::string &s) {
        if (l <= level) {
            sprint_nolr(ostream, s);
        }
    }

    template <Level l>
    inline void smprint(std::ostream &ostream, const std::string &s) {
        if (l <= level && my_rank == io_rank) {
            sprint_nolr(ostream, s);
        }
    }

    template <Level l>
    inline void scprint(std::ostream &ostream, const std::string &s) {
        if (l <= level) {
            scprint_nol(ostream, s);
        }
    }

  protected:
    IOMpi();

    ~IOMpi();

    int Get_io_rank();

    void scprint_nol(std::ostream &ostream, const std::string &s);

    void sprint_nolr(std::ostream &ostream, const std::string &s);

    static const int DEFAULT_IO_PROC = 0;

    void *void_io_comm;
    int io_rank, my_rank, comm_size;
    int cprintf_tag;
    pthread_mutex_t cprintf_mutex;
    Level level;
};

#endif