# hp-Adaptive Signed Distance Fields
A multi-threaded implementation of [An hp-Adaptive Discretization Algorithm for
Signed Distance Field Generation](https://www.animation.rwth-aachen.de/media/papers/2017-TVCG-HPDistanceFields.pdf). Through an octree, the library generates an approximation of any function (and simple operations between such functions) defined over an arbitrary domain. Parameters controlling quality, thread usage and whether to perform a continuity post-process are all fully exposed to the user. At the cost of a small approximation error, query perfomance for relatively complex input SDFs can be improved by orders of magnitude. 

<figure>
    <img src="https://i.imgur.com/HuNFjhh.jpeg"
         alt="Some SDFs">
    <figcaption>SDF approximations with varying target errors and pre/post continuity optimisation. Note the relatively small visual improvement on the top row as the target error is very small.</figcaption>
</figure>

### Library API

All library objects and functions are namespaced within SDF. The API is intentionally very simple and creating an approximation can be done in a dozen or so lines of code, as shown below:

```
  // Example SDF
  auto SphereFunc = [](const Eigen::Vector3d& pt_) -> double
  {
      return pt_.norm() - 0.25;
  };
  
  // Create config
  SDF::Config config;
  config.targetErrorThreshold       = 0.000001;
  config.nearnessWeighting.type     = SDF::Config::NearnessWeighting::Type::Exponential;
  config.nearnessWeighting.strength = 3.0;
  config.continuity.enforce         = true;
  config.continuity.strength        = 8.0;
  config.threadCount                = 12;
  config.root                       = Eigen::AlignedBox3f(Eigen::Vector3f(-0.25, -0.25, -0.25), Eigen::Vector3f(0.5, 0.5, 0.5));
  
  // Create tree and pass config and SDF
  SDF::Octree hpOctree;
  hpOctree.Create(config, SphereFunc)
  
  // Query the tree as required
  const Eigen::Vector3d somePt = Eigen::Vector3d::Zero();
  const double someVal         = hpOctree.Query(somePt);
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

### Examples

For each example, we have a slice of the approximated SDF and the octree structure at z = 0. In these cases, the function being approximated was a simple union for shapes that have closed-form SDFs readily available on the Internet. The colours in the right portion of the image correspond to basis polynomial degree with:

* Grey        = 2
* Green       = 3
* Light Blue  = 4
* Medium Blue = 5
* Dark Blue   = 6

To add context to the perfomance figures, all examples were generated on a PC equipped with a Ryzen 5800X and OpenMP enabled to speed up Eigen's ConjugateGradient solver. 

![ImageA](https://i.imgur.com/HPIm2IM.png)

![ImageB](https://i.imgur.com/W1MgeER.png)

![ImageC](https://i.imgur.com/ZW1iW9y.png)

### Future Improvements

* Bite the bullet and significantly improve the CMake code and test compilation on other platforms (Linux etc).
* Bug somewhere in the Meshing code with the PseudoNormal calculation. That is, Mesh.FXYZ(pt) can sometimes return the wrong sign and will cause the SDF approximation to fail or look a bit silly. This is quite rare (generated this video from a 2Mil trimesh https://www.youtube.com/watch?v=8LauGgRzlcI) but worth considering. 
