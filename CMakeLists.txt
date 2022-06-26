cmake_minimum_required(VERSION 3.0)

project(hpAdaptiveSDF VERSION 1.0)

# C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# OpenMP for Eigen
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Code includes
include_directories("External/eigen-master")
include_directories("External/stb-master")
add_subdirectory("Source/HP")
add_subdirectory("Source/Meshing")

# Add exe
add_executable(Example "Source/Main.cpp")