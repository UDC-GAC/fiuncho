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
 * @file IOMpi.h
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief IOMpi class definition.
 */

#ifndef MPI3SNP_IOMPI_H
#define MPI3SNP_IOMPI_H

#include <cstdarg>
#include <iostream>
#include <pthread.h>

/**
 * @class IOMpi
 *
 * @brief Handles all the I/O throughout the program. Uses Meyers' Singleton pattern, which describes a thread-safe
 * singleton initialization (C++11 or superior). The object destruction is handled at the end of the program.
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

    static void Set_print_level(Level l) {
        Instance().level = l;
    }

    // delete copy and move constructors and assign operators
    IOMpi(IOMpi const &) = delete;             // Copy construct
    IOMpi(IOMpi &&) = delete;                  // Move construct
    IOMpi &operator=(IOMpi const &) = delete;  // Copy assign
    IOMpi &operator=(IOMpi &&) = delete;      // Move assign

    template<Level l>
    inline void print(const std::string &s) {
        if (l <= level) {
            sprint_nolr(std::cout, s);
        }
    }

    template<Level l>
    inline void mprint(const std::string &s) {
        if (l <= level && my_rank == io_rank) {
            sprint_nolr(std::cout, s);
        }
    }

    template<Level l>
    inline void cprint(const std::string &s) {
        if (l <= level) {
            scprint_nol(std::cout, s);
        }
    }

    template<Level l>
    inline void sprint(std::ostream &ostream, const std::string &s) {
        if (l <= level) {
            sprint_nolr(ostream, s);
        }
    }

    template<Level l>
    inline void smprint(std::ostream &ostream, const std::string &s) {
        if (l <= level && my_rank == io_rank) {
            sprint_nolr(ostream, s);
        }
    }

    template<Level l>
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

#endif //MPI3SNP_IOMPI_H