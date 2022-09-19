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