#pragma once

#include <vector>
#include <queue>
#include <functional>

#include <malloc.h>

#include "Eigen/Core"
#include "Eigen/Geometry"
#include "Eigen/Sparse"

#define HAS_STB 0
#if HAS_STB
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#endif

#include "Literals.h"

#include "Utility.h"
#include "Legendre.h"
#include "Node.h"
#include "Config.h"
#include "BuildThreadPool.h"
#include "ContinuityThreadPool.h"
#include "MemoryBlock.h"

namespace SDF
{
	class Octree
	{
	public:
		Octree();
		~Octree();

		/// Don't want to allow copying just yet
		Octree& operator=(const Octree& other_) = delete;
		Octree& operator=(Octree&& other_) = delete;

		/// Approximates F_ inside [-0.5f, 0.5f]^3 using the parameters in config_
		void Create(const Config& config_, std::function<f64(const Eigen::Vector3d& pt_)> F_);

		/// Resets the tree
		void Clear();

		/// Creates an octree from a previously serialised version
		void FromMemoryBlock(MemoryBlock octBlock_);

		/// Serialises an octree to a memory block owned by malloc
		MemoryBlock ToMemoryBlock();

		/// Returns the approximated distance from F = 0. Negative implies inside the boundary
		f64 Query(const Eigen::Vector3d& pt_);

		/// Writes a 2048*2048 image of the approximated SDF about z = c_ to fName_.bmp
#if HAS_STB
		void OutputFunctionSlice(const char* fName_, const f64 c_);
#endif

	private:
		friend class BuildThreadPool;
		friend class ContinuityThreadPool;

		/// Predicate for use in the priority queue
		struct PriorityQueuePredicate
		{
			bool operator()(const std::pair<usize, f64>& v1_, const std::pair<usize, f64>& v2_)
			{
				return v1_.second < v2_.second;
			}
		};
		std::priority_queue<std::pair<usize, f64>, std::vector<std::pair<usize, f64>>, PriorityQueuePredicate> nodeQueue;

		std::map<std::pair<usize, usize>, bool> procMap;
		Config config;
		std::function<f64(const Eigen::Vector3d& pt_)> F;
		f64* coeffStore;
		std::vector<Node> nodes;

		/// Estimates the improvement we'd obtain by subviding the node
		f64 EstimateHImprovement(const BuildThreadPool::Input& inputData_, Node::Basis* hBases_, f64* hBasesErrors_);

		/// Estimates the improvement we'd obtain by increasing the node's basis degree
		f64 EstimatePImprovement(const BuildThreadPool::Input& inputData_, Node::Basis& pBasis_, f64& pBasisError_);

		/// Fits a basis of degree_ to basis_ over the volume aabb_. Returns the error
		f64 FitPolynomial(Node::Basis& basis_, const Eigen::AlignedBox3f& aabb_, const usize degree_, const usize depth_);

		/// Subdivides nodes[nodeIdx_]. Doesn't fit any polynomials to the children
		void Subdivide(const usize nodeIdx_);

		/// Returns the ith corner AABB of aabb_
		Eigen::AlignedBox3f CornerAABB(const Eigen::AlignedBox3f& aabb_, const usize i_);

		/// Creates a subdivided root node
		void CreateRoot();

		/// Manages the threadpool responsible for building the tree
		void RunBuildThreadPool(f64& totalCoeffError_);

		/// Evaluates (4) for a given p_ and x_
		f64 LpX(const usize p_, const f64 x_);

		/// Evaluates (2) using the given basis_ at pt_
		f64 FApprox(const Node::Basis& basis_, const Eigen::AlignedBox3f& aabb_, const Eigen::Vector3d& pt_, const usize depth_);

		/// Applies (11)
		f64 CalculatePolyWeighting(const Node::Basis& basis_, const Eigen::AlignedBox3f& aabb_, const usize depth_);

		/// Applies (12)
		f64 CalculateExpWeighting(const Node::Basis& basis_, const Eigen::AlignedBox3f& aabb_, const usize depth_);

		/// Uniformly divides the octree and fits a polynomial of degree 2 at every depth 4 node
		void UniformlyRefine(f64& totalCoeffError_);

		/// Performs the main algo loop
		void TickBuildThread(const BuildThreadPool::InitialData* const threadData_);

		/// Performs the continuity main algo loop
		void TickContinuityThread(const ContinuityThreadPool::InitialData* const threadData_, const usize threadIdx_);

		/// Performs a continuity post process to the SDF to try and get rid of discontinuities between cells
		void PerformContinuityPostProcess(const usize nCoeffs_);

		/// Cleans up the memory by making it s.t. each node's coeff ptr is an offset into the same large block of memory. ~10% speed increase
		usize ReallocCoeffs();

		/// Proc functions from https://www.cs.rice.edu/~jwarren/papers/dualcontour.pdf
		void NodeProc(const usize nodeIdx_, std::queue<ContinuityThreadPool::Input>& jobQueue_);
		void FaceProcX(const usize nodeIdxA_, const usize nodeIdxB_, std::queue<ContinuityThreadPool::Input>& jobQueue_);
		void FaceProcY(const usize nodeIdxA_, const usize nodeIdxB_, std::queue<ContinuityThreadPool::Input>& jobQueue_);
		void FaceProcZ(const usize nodeIdxA_, const usize nodeIdxB_, std::queue<ContinuityThreadPool::Input>& jobQueue_);

		/// Fills a triplet array for the sparse matrix up via another threadpool (if nThreads > 1)
		void RunContinuityThreadPool(std::vector<Eigen::Triplet<f64>>& integralMatrixTriplets_);

		/// Evaluates the shared face integral as in Appendix A... Though differs slightly due to possible error in paper
		void EvaluateSharedFaceIntegralAnalytically(const usize nodeIdxA_, const usize nodeIdxB_, const u8 dim_, std::vector<Eigen::Triplet<f64>>& matTriplets_);

		/// Evaluates the shared face integral using standard GL quadrature
		void EvaluateSharedFaceIntegralNumerically(const usize nodeIdxA_, const usize nodeIdxB_, const u8 dim_, std::vector<Eigen::Triplet<f64>>& matTriplets_);
	};
}