# CMake setup
cmake_minimum_required(VERSION 3.11)

# Project setup
project(hpAdaptiveSDF VERSION 1.1)

# C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories("Include")
include_directories("External/eigen")
include_directories("External/stb")

# Libs and exe
add_subdirectory(Source)
add_subdirectory(Source/Tests)