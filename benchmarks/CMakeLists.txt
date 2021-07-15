# This file is part of Fiuncho.
#
# Fiuncho is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Fiuncho is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Fiuncho. If not, see <http://www.gnu.org/licenses/>.

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
