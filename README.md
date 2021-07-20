# hp-Adaptive Signed Distance Fields

A multithreaded implementation of [An hp-Adaptive Discretization Algorithm for
Signed Distance Field Generation](https://www.animation.rwth-aachen.de/media/papers/2017-TVCG-HPDistanceFields.pdf). The code generates an approximation of an abitrary SDF defined over the unit cube. A example use case is to quickly determine the signed distance of a point from a triangular mesh, which can be quite costly even when using a BVH. At the cost of small approximation error, query perfomance can be improved by orders of magnitude. An example use of the API is given below:

```
  // Example SDF
  auto SphereFunc = [](const Eigen::Vector3d& pt_) -> double
  {
      return pt_.norm() - 0.25;
  };
  
  // Create config
  SDF::Config config;
  config.maximumSDFError = 0.000001;
  config.nearnessWeighting.type = SDF::Config::NearnessWeighting::Type::Exponential;
  config.nearnessWeighting.stength = 3.0;
  config.continuity.enforce = true;
  config.continuity.strength = 8.0;
  config.threadCount = 12;
  
  // Create tree and pass config and SDF
  SDF::Octree hpOctree;
  hpOctree.Create(config, SphereFunc)
  
  // Query the tree as required
  const Eigen::Vector3d somePt = Eigen::Vector3d::Zero();
  const double someVal = hpOctree.Query(somePt);
```

For more complex SDFs and smaller target errors, serialisation is fully supported via SDF::MemoryBlock, which is a simple pointer/size struct owned by malloc:

```
  SDF::MemoryBlock hpOctreeBlock = hpOctree.ToMemoryBlock();
  
  /*
      ...
      Move block around memory or write to disk
      ...
  */
  
  // Free block pointer
  free(hpOctreeBlock.ptr);
  
  /*
      ...
      Obtain pointer and size of block from disk or memory
      ...
  */
  
  SDF::Octree someNewOctree;
  someNewOctree.FromMemoryBlock(SDF::MemoryBlock{ someSize, somePtr });
```

Some octree examples are given below. For each one, we have a slice of the approximated SDF and the octree structure at z = 0. One thing to note is the effect of the continuity optimisation between images B and C, which otherwise share the exact same config. The colours in the right portion of the image correspond to basis polynomial degree with:

* Grey = 2
* Green = 3
* Light Blue = 4
* Medium Blue = 5
* Dark Blue = 6

All examples were generated using a Ryzen 5800X and OpenMP enabled to speed up Eigen's ConjugateGradient solver.

![ImageA](https://i.imgur.com/kCpQJvk.png)

![ImageB](https://i.imgur.com/cbhtWNn.png)

![ImageC](https://i.imgur.com/0aHdBCW.png)

An ordered ToDo list can be found in ToDo.txt (whatever else).
