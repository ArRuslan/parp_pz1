cmake_minimum_required(VERSION 3.29)
project(parp_pz1)

set(CMAKE_CXX_STANDARD 17)

find_package(OpenMP REQUIRED)

add_executable(parp_pz1 main.cpp)
target_link_libraries(parp_pz1 PRIVATE OpenMP::OpenMP_CXX)
target_link_options(parp_pz1 PRIVATE -Wl,-z,stack-size=67108864)