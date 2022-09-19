#pragma once

#include <vector>
#include <queue>
#include <functional>
#include <map>
#include <limits>
#include <mutex>
#include <atomic>
#include <thread>

#include <malloc.h>

#define EIGEN_MPL2_ONLY 
#include "Eigen/Core"
#include "Eigen/Geometry"
#include "Eigen/Sparse"

#include "Utility/Literals.h"
#include "Utility/MemoryBlock.h"
#include "HP/Consts.h"
#include "HP/Ray.h"
#include "HP/Utility.h"
#include "HP/Legendre.h"
#include "HP/Node.h"
#include "HP/Config.h"
#include "HP/BuildThreadPool.h"
#include "HP/ContinuityThreadPool.h"

namespace SDF
{
	class Octree
	{
	public:
		Octree();
		~Octree();

		/// Copying/Moving
        Octree(const Octree& other_);
        Octree(Octree&& other_);
		Octree& operator=(const Octree& other_);
		Octree& operator=(Octree&& other_);

		/// Approximates F_ using the parameters in config_
		void Create(const Config& config_, std::function<f64(const Eigen::Vector3d& pt_, const u32 threadIdx_)> F_);

        /// Resultant SDF = Min(oldF, F_)
        void UnionSDF(std::function<f64(const Eigen::Vector3d& pt_, const u32 threadIdx_)> F_);

        /// Resultant SDF = Max(-oldF, F_)
        void SubtractSDF(std::function<f64(const Eigen::Vector3d& pt_, const u32 threadIdx_)> F_);

        /// Resultant SDF = Max(oldF, F_)
        void IntersectSDF(std::function<f64(const Eigen::Vector3d& pt_, const u32 threadIdx_)> F_);

		/// Resets the tree
		void Clear();

		/// Creates an octree from a previously serialised version
		void FromMemoryBlock(MemoryBlock octBlock_);

		/// Serialises an octree to a memory block owned by malloc
		MemoryBlock ToMemoryBlock() const;

		/// Returns the approximated distance from F = 0. Negative implies inside the boundary
		f64 Query(const Eigen::Vector3d& pt_) const;

        /// NOTE: Untested
        /// Returns whether a ray intersected the SDF and the value of the first hit if it does
        bool QueryRay(const Ray& ray_, const f64 tMax_, f64& t_) const;

		/// As with Query, but with an optional unit gradient calculated via CD
		f64 QueryWithGradient(const Eigen::Vector3d& pt_, Eigen::Vector3d& unitNormal_) const;

        /// Returns the aabb of the root node
        Eigen::AlignedBox3f GetRootAABB() const;

#ifdef HAS_STB
		/// Writes a 2048*2048 image of the approximated SDF about z = c_ to fName_.bmp
		void OutputFunctionSlice(const char* fName_, const f64 c_, const Eigen::AlignedBox3f& viewArea_);
#endif

	private:
        static constexpr f64 INITIAL_NODE_ERR = 100.0;

		friend class BuildThreadPool;
		friend class ContinuityThreadPool;

		/// Predicate for use in the priority queue
		struct PriorityQueuePredicate
		{
			bool operator()(const std::pair<u32, f64>& v1_, const std::pair<u32, f64>& v2_)
			{
				return v1_.second < v2_.second;
			}
		};
		std::priority_queue<std::pair<u32, f64>, std::vector<std::pair<u32, f64>>, PriorityQueuePredicate> nodeQueue;

		std::map<std::pair<u32, u32>, bool>                                  procMap;
		std::function<f64(const Eigen::Vector3d& pt_, const u32 threadIdx_)> F;
		f64*                                                                 coeffStore;
		std::vector<Node>                                                    nodes;

        Config          config;
        Eigen::Vector3d configRootCentre;
        Eigen::Vector3d configRootInvSizes;

		/// Estimates the improvement we'd obtain by subviding the node
		f64 EstimateHImprovement(const BuildThreadPool::Input& inputData_, Node::Basis* hBases_, f64* hBasesErrors_, const u32 threadIdx_);

		/// Estimates the improvement we'd obtain by increasing the node's basis degree
		f64 EstimatePImprovement(const BuildThreadPool::Input& inputData_, Node::Basis& pBasis_, f64& pBasisError_, const u32 threadIdx_);

		/// Fits a basis of degree_ to basis_ over the volume aabb_. Returns the error
		f64 FitPolynomial(Node::Basis& basis_, const Eigen::AlignedBox3f& aabb_, const u8 degree_, const u32 depth_, const u32 threadIdx_);

		/// Subdivides nodes[nodeIdx_]. Doesn't fit any polynomials to the children
		void Subdivide(const u32 nodeIdx_);

		/// Returns the ith corner AABB of aabb_
		Eigen::AlignedBox3f CornerAABB(const Eigen::AlignedBox3f& aabb_, const u32 i_);

		/// Creates a subdivided root node
		void CreateRoot();

		/// Manages the threadpool responsible for building the tree
		void RunBuildThreadPool();

		/// Evaluates (4) for a given p_ and x_
		f64 LpX(const u32 p_, const f64 x_);

		/// Evaluates (2) using the given basis_ at pt_
		f64 FApprox(const Node::Basis& basis_, const Eigen::AlignedBox3f& aabb_, const Eigen::Vector3d& pt_, const u32 depth_) const;
		
		/// Same with FApprox but with optimisations to allow efficient gradient computation
		f64 FApproxWithGradient(const Node::Basis& basis_, const Eigen::AlignedBox3f& aabb_, const Eigen::Vector3d& pt_, const u32 depth_, Eigen::Vector3d& unitNormal_) const;

		/// Applies (11)
		f64 CalculatePolyWeighting(const Node::Basis& basis_, const Eigen::AlignedBox3f& aabb_, const u32 depth_);

		/// Applies (12)
		f64 CalculateExpWeighting(const Node::Basis& basis_, const Eigen::AlignedBox3f& aabb_, const u32 depth_);

		/// Uniformly divides the octree and fits a polynomial of degree 2 at every depth 4 node
		void UniformlyRefine();

		/// Performs the main algo loop
		void TickBuildThread(const BuildThreadPool::InitialData* const threadData_, const u32 threadIdx_);

		/// Performs the continuity main algo loop
		void TickContinuityThread(const ContinuityThreadPool::InitialData* const threadData_, const u32 threadIdx_);

		/// Performs a continuity post process to the SDF to try and get rid of discontinuities between cells
		void PerformContinuityPostProcess(const u32 nCoeffs_);

		/// Cleans up the memory by making it s.t. each node's coeff ptr is an offset into the same large block of memory. ~10% speed increase
		u32 ReallocCoeffs();

		/// Proc functions from https://www.cs.rice.edu/~jwarren/papers/dualcontour.pdf
		void NodeProc(const u32 nodeIdx_, std::queue<ContinuityThreadPool::Input>& jobQueue_);
		void FaceProc(const u32 nodeIdxA_, const u32 nodeIdxB_, const u8 dim_, std::queue<ContinuityThreadPool::Input>& jobQueue_);

		/// Fills a triplet array for the sparse matrix up via another threadpool (if nThreads > 1)
		void RunContinuityThreadPool(std::vector<Eigen::Triplet<f64>>& integralMatrixTriplets_);

		/// Evaluates the shared face integral as in Appendix A... Though differs slightly due to possible error in paper
		void EvaluateSharedFaceIntegralAnalytically(const u32 nodeIdxA_, const u32 nodeIdxB_, const u8 dim_, std::vector<Eigen::Triplet<f64>>& matTriplets_);

		/// Evaluates the shared face integral using standard GL quadrature
		void EvaluateSharedFaceIntegralNumerically(const u32 nodeIdxA_, const u32 nodeIdxB_, const u8 dim_, std::vector<Eigen::Triplet<f64>>& matTriplets_);
	};
}