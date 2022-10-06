#pragma once

#include <queue>
#include <mutex>
#include <atomic>
#include <thread>

#include <malloc.h>

#include "Eigen/Core"
#include "Eigen/Geometry"
#include "Eigen/Sparse"

#include "Utility/Literals.h"
#include "HP/Config.h"

namespace SDF
{
    class Octree;

    class ContinuityThreadPool
    {
    public:
        struct Input
        {
            std::pair<u32, u32> nodeIdxs;
            u8                  dim;
        };

        struct InitialData
        {
            Config config;

            std::atomic<bool>* shutdownAtom;

            std::queue<Input>* inputQueue;
            std::mutex*        inputQueueMutex;

            std::vector<Eigen::Triplet<f64>>* matTriplets;
        };

        void StartThreads(InitialData* const initialData_, Octree* const octree_);
        void StopThreads();

    private:
        std::vector<std::thread> threadPool;
    };
}
