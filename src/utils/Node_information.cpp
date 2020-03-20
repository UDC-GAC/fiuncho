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
 * @file Gpu_node_information.cpp
 * @author Christian Ponte
 * @date 05 Feb 2018
 *
 * @brief Implementation of abstract class Node_information members to be
 * inherited.
 */

#include <algorithm>
#include <array>
#include <fiuncho/Definitions.h>
#include <fiuncho/engine/cpu/Cpu_node_information.h>
#include <fiuncho/utils/Node_information.h>
#include <fstream>
#include <mpi.h>

Node_information *Node_information::Builder::get_information() {
    return new Cpu_node_information();
}

Node_information *Node_information::build_from_byteblock(const void *ptr) {
    return new Cpu_node_information(ptr);
}

Node_information::Node_information() {
    std::ifstream ifs(hardware_id_file);
    std::string content((std::istreambuf_iterator<char>(ifs)),
                        (std::istreambuf_iterator<char>()));
    std::hash<std::string> hash_function;
    hardware_id = std::to_string(hash_function(content));
}

std::vector<Node_information *> Node_information::gather(int process) {
    int rank, count;

    MPI_Comm_size(MPI_COMM_WORLD, &count);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Node_information *local = Node_information::Builder::get_information();
    std::vector<Node_information *> output;

    if (rank == process) {
        output.push_back(local);

        MPI_Status status;
        char *block;
        int block_size;
        Node_information *temp;
        for (int i = 0; i < count; i++) {
            if (i != process) {
                MPI_Probe(i, 0, MPI_COMM_WORLD, &status);
                MPI_Get_count(&status, MPI_BYTE, &block_size);
                block = new char[block_size];
                MPI_Recv(block, block_size, MPI_BYTE, i, 0, MPI_COMM_WORLD,
                         nullptr);
                temp = Node_information::build_from_byteblock(block);
                delete[] block;
                auto pos = std::find_if(output.begin(), output.end(),
                                        [&temp](Node_information *it) {
                                            return it->hardware_id ==
                                                   temp->hardware_id;
                                        });
                if (pos == output.end()) {
                    output.push_back(temp);
                } else {
                    (*pos)->add_processes(temp->processes());
                    delete temp;
                }
            }
        }
    } else {
        char *block;
        int size = local->to_byteblock((void **)&block);
        MPI_Send(block, size, MPI_BYTE, process, 0, MPI_COMM_WORLD);
        delete[] block;
        delete local;
    }
    return output;
}
