# CMake setup
cmake_minimum_required(VERSION 3.11)

# Project setup
project(hpAdaptiveSDF VERSION 1.0)

# C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# OpenMP for Eigen
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories("Include")
include_directories("External/eigen-master")
include_directories("External/stb-master")

# .obj resources
file(COPY Resources DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Add exe
add_subdirectory("Source")
add_executable(Tests "Source/Main.cpp")