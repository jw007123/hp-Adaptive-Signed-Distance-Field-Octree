#ifdef _WIN32
#define _CRT_SECURE_NO_WARNINGS 1
#endif

#include <stdio.h>

#include "HPUnitTests.cpp"
#include "HPBenchmarks.cpp"
#include "MeshingUnitTests.cpp"
#include "MeshingBenchmarks.cpp"

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