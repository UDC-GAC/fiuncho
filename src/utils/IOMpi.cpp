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
 * @file IOMpi.cpp
 * @author Christian Ponte
 * @date 1 March 2018
 *
 * @brief IOMpi class members implementation.
 */

#include "IOMpi.h"
#include <mpi.h>

// Macro to cast void_io_comm attribute from type void * into MPI_Comm type for better readability
#ifndef io_comm
#define io_comm *(MPI_Comm *) void_io_comm
#endif

IOMpi::IOMpi() {
    void_io_comm = new MPI_Comm;
    MPI_Comm_dup(MPI_COMM_WORLD, &io_comm);
    MPI_Comm_size(io_comm, &comm_size);
    MPI_Comm_rank(io_comm, &my_rank);
    io_rank = Get_io_rank();
    cprintf_tag = 0;
    pthread_mutex_init(&cprintf_mutex, nullptr);
    level = N;
}

IOMpi::~IOMpi() {
    int finalized;
    MPI_Finalized(&finalized);
    if (!finalized)
        MPI_Comm_free((MPI_Comm *) io_comm);
    pthread_mutex_destroy(&cprintf_mutex);
}

int IOMpi::Get_io_rank() {
    // Determine IO process
    int mpi_io, flag, io_rank;

    MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_IO, &mpi_io, &flag);
    if (!flag) {
        // Attribute not cached
        io_rank = DEFAULT_IO_PROC;
    } else if (mpi_io == MPI_PROC_NULL || mpi_io == MPI_ANY_SOURCE || comm_size < mpi_io) {
        io_rank = DEFAULT_IO_PROC;
    } else {
        // Multiple IO processes
        int err = MPI_Allreduce(&mpi_io, &io_rank, 1, MPI_INT, MPI_MIN, io_comm);
        if (err != MPI_SUCCESS) {
            int error_class, len;
            char error_string[1024];
            MPI_Error_class(err, &error_class);
            MPI_Error_string(error_class, error_string, &len);
            std::cerr << error_string << std::endl;
            io_rank = DEFAULT_IO_PROC;
        } else if (io_rank < 0 || io_rank >= comm_size) {
            io_rank = DEFAULT_IO_PROC;
        }
    }
    return io_rank;
}


void IOMpi::scprint_nol(std::ostream &ostream, const std::string &s) {
    int tag;

    pthread_mutex_lock(&cprintf_mutex);
    tag = cprintf_tag++;
    pthread_mutex_unlock(&cprintf_mutex);

    if (my_rank == io_rank) {
        MPI_Message message;
        MPI_Status status;
        int count;
        for (int i = 0; i < comm_size; i++) {
            if (i == my_rank) {
                ostream << "Process " + std::to_string(my_rank) + " > " << s;
            } else {
                if (MPI_Mprobe(i, tag, io_comm, &message, &status) != MPI_SUCCESS) {
                    return;
                }
                MPI_Get_count(&status, MPI_CHAR, &count);
                char str[count];
                if (MPI_Mrecv(str, count, MPI_CHAR, &message, nullptr) != MPI_SUCCESS) {
                    return;
                }
                ostream << "Process " + std::to_string(i) + " > " << str;
            }
        }
        std::flush(ostream);
    } else {
        if (MPI_Send(s.c_str(), s.length() + 1, MPI_CHAR, io_rank, tag, io_comm) != MPI_SUCCESS) {
            return;
        }
    }
}

void IOMpi::sprint_nolr(std::ostream &ostream, const std::string &s) {
    ostream << s;
    std::flush(ostream);
}