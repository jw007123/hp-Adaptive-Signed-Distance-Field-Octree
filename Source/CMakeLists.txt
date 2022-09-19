# OpenMP for Eigen
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_library(Lib_HP HP/BuildThreadPool.cpp
				   HP/Config.cpp
				   HP/ContinuityThreadPool.cpp
				   HP/Node.cpp
				   HP/Octree.cpp
				   HP/Ray.cpp
				   Meshing/Utility.cpp
				   Meshing/BVH.cpp
				   Meshing/Mesh.cpp
				   Meshing/NNOctree.cpp
				   Meshing/ObjParser.cpp) 

target_include_directories (Lib_HP PUBLIC "../Include")
target_include_directories (Lib_HP PUBLIC "../External/eigen-master")
target_include_directories (Lib_HP PUBLIC "../External/stb-master")