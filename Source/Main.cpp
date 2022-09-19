#ifdef _WIN32
#define _CRT_SECURE_NO_WARNINGS 1
#endif

#include <stdio.h>

#include "HP/Ray.cpp"
#include "HP/Config.cpp"
#include "HP/Node.cpp"
#include "HP/Octree.cpp"
#include "HP/ContinuityThreadPool.cpp"
#include "HP/BuildThreadPool.cpp"
#include "HP/UnitTests.cpp"
#include "HP/Benchmarks.cpp"

#include "Meshing/Utility.cpp"
#include "Meshing/ObjParser.cpp"
#include "Meshing/NNOctree.cpp"
#include "Meshing/BVH.cpp"
#include "Meshing/Mesh.cpp"
#include "Meshing/UnitTests.cpp"
#include "Meshing/Benchmarks.cpp"

bool RunTests()
{
    // HP
    {
        printf("Running HP unit tests...\n\n");

        SDF::UnitTests unitTests;
        if (!unitTests.Run())
        {
            std::this_thread::sleep_for(std::chrono::seconds(5));
            return false;
        }
    }

    // Meshing
    {
        printf("Running Meshing unit tests...\n\n");

        Meshing::UnitTests unitTests;
        if (!unitTests.Run())
        {
            std::this_thread::sleep_for(std::chrono::seconds(5));
            return false;
        }
    }
    
    std::this_thread::sleep_for(std::chrono::seconds(5));
    return true;
}
   

void RunBenchmarks()
{
    // HP
    {
        printf("Running HP benchmarks...\n\n");

        SDF::Benchmarks benchmarks;
        benchmarks.Run();

        std::this_thread::sleep_for(std::chrono::seconds(5));
    }

    // Meshing
    {
        printf("Running Meshing benchmarks...\n\n");

        Meshing::Benchmarks benchmarks;
        benchmarks.Run();

        std::this_thread::sleep_for(std::chrono::seconds(5));
    }
}


int main()
{
    if (RunTests())
    {
        RunBenchmarks();
    }

	return 0;
}