#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>

#include "HP/Ray.cpp"
#include "HP/Config.cpp"
#include "HP/Node.cpp"
#include "HP/Octree.cpp"
#include "HP/ContinuityThreadPool.cpp"
#include "HP/BuildThreadPool.cpp"
#include "HP/UnitTests.cpp"
#include "HP/Benchmarks.cpp"

#include "Meshing/Simplex.cpp"
#include "Meshing/ObjParser.cpp"
#include "Meshing/NNOctree.cpp"
#include "Meshing/BVH.cpp"
#include "Meshing/Mesh.cpp"
#include "Meshing/Benchmarks.cpp"

void RunTests()
{
    // HP
    {
        printf("Running HP unit tests...\n\n");

        SDF::UnitTests unitTests;
        if (!unitTests.Run())
        {
            std::this_thread::sleep_for(std::chrono::seconds(5));
            return;
        }
    }
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


void Create()
{
    Meshing::Mesh objMesh;
    objMesh.CreateFromObj("C:\\VS2017Projects\\SDF\\Mesh-BVH-for-hp\\Resources\\crank.obj");
    if (!objMesh.HasValidState())
    {
        return;
    }

    Meshing::BVH objBVH;
    objBVH.Create(objMesh);
    if (!objBVH.HasValidState())
    {
        return;
    }

    Eigen::AlignedBox3f rootAABB = objMesh.CalculateMeshAABB();
    {
        f32 maxEdge = 0.0;
        for (u8 i = 0; i < 3; ++i)
        {
            maxEdge = std::max<f32>(maxEdge, (rootAABB.max()(i) - rootAABB.min()(i)));
        }
        for (u8 i = 0; i < 3; ++i)
        {
            rootAABB.max()(i) = rootAABB.min()(i) + maxEdge * 1.1;
            rootAABB.min()(i) = rootAABB.min()(i) - maxEdge * 0.1;
        }
    }

    auto MeshFunc = [&objBVH, &objMesh](const Eigen::Vector3d& pt_, const u32 threadIdx_) -> f64
    {
        return objMesh.SignedDistanceAtPt(pt_.cast<f32>(), objBVH, threadIdx_);
    };

    SDF::Config hpConfig;
    hpConfig.targetErrorThreshold = pow(10, -7);
    hpConfig.nearnessWeighting.type = SDF::Config::NearnessWeighting::Exponential;
    hpConfig.nearnessWeighting.strength = 3.0;
    hpConfig.continuity.enforce = true;
    hpConfig.continuity.strength = 8.0;
    hpConfig.threadCount = 16;
    hpConfig.root = rootAABB;
    hpConfig.enableLogging = true;

    SDF::Octree hpOctree;
    hpOctree.Create(hpConfig, MeshFunc);

    MemoryBlock blk = hpOctree.ToMemoryBlock();
    HANDLE hFile = CreateFile("Crank.bin", GENERIC_WRITE, FILE_SHARE_READ, NULL, CREATE_NEW, FILE_ATTRIBUTE_NORMAL, NULL);
    if (hFile != INVALID_HANDLE_VALUE)
    {
        DWORD bytesWritten;
        WriteFile(hFile, blk.ptr, blk.size, &bytesWritten, nullptr);
        CloseHandle(hFile);
    }
    free(blk.ptr);

    Eigen::AlignedBox3f r = hpOctree.GetRootAABB();
    r.min() += Eigen::Vector3f(0.1f, 0.1f, 0.1f);
    r.max() -= Eigen::Vector3f(0.1f, 0.1f, 0.1f);
    u32 cnt = 0;
    for (f32 i = -0.5f; i < 0.5f; i += 0.0025f)
    {
        char s[128];
        sprintf(s, "slice-%zu", cnt);
        hpOctree.OutputFunctionSlice(s, i, r);
        cnt++;
    }
}


void Read()
{
    HANDLE hFile = CreateFile("Ramesses.bin", GENERIC_READ, 0, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);

    MemoryBlock blk;
    blk.ptr = malloc(1024 * 1024 * 6);
    blk.size = 1024 * 1024 * 6;
    DWORD bytesRead;
    ReadFile(hFile, blk.ptr, blk.size, &bytesRead, nullptr);
    CloseHandle(hFile);
    blk.ptr = realloc(blk.ptr, bytesRead);
    blk.size = bytesRead;

    SDF::Octree hpOctree;
    hpOctree.FromMemoryBlock(blk);
    free(blk.ptr);

    Eigen::AlignedBox3f r = hpOctree.GetRootAABB();
    r.min() += Eigen::Vector3f(0.1f, 0.1f, 0.1f);
    r.max() -= Eigen::Vector3f(0.1f, 0.1f, 0.1f);
    u32 cnt = 0;
    for (f32 i = -0.5f; i < 0.5f; i += 0.0025f)
    {
        char s[128];
        sprintf(s, "slice-%zu", cnt);
        hpOctree.OutputFunctionSlice(s, i, r);
        cnt++;
    }
}


int main()
{
    Create();

	return 0;
}