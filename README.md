#**landing-optimization-framework**

This file generates a main() to optimize the control parameters (desired length
length gain, speed gain) of muscle stretch controllers to match a given motion.
This will build an executable that can take the following arguments:

*  -m <optimizedDropModel.osim>: model file to be loaded
*  -csfile <optimizedControllerSet.xml>: initial guess of controller set parameters for optimizer
*  -kfile <DesiredKinematics.sto>: trajectories of desired kinematics
*  -ti <0.0> default time from which to start simulating
*  -tf <0.1> default simulation end time
*  -opt: flag to run optimization
*  -resume: flag to run optimization but resume using resumecmaes.dat file
*  -lambda <int>: CMA parameter, set number of samples per outer iteration
*  -sigma <double>: CMA parameter, set initial constant of covariance matrix
*  -maxIters <int>: optimizer parameter, set max outer iterations

###Usage: 
```
dropLandingSim_nodal.exe -m InputDropModel.osim -csfile InitialControllerSet.xml -kfile DesiredKinematics.sto -ti 0.0 -tf 0.1 -opt -lambda 50 -sigma 0.1 -maxIters 2000
```
The -opt flag is used when an optimization should be run. If the -opt flag is
not present, then a single forward simulation with analyses will be performed.

Parallelization of the CMA optimizer is done by using MPI. If ENABLE_MPI was
checked in CMake, then you can call on the executable to use more threads using:
```
mpiexec -n <numthreads> dropLandingSim_nodal.exe -opt -lambda 50 ...
```

###Dependencies
1. OpenSim 3.2 or above by [installing a distribution](https://simtk.org/home/opensim) or [building from source](https://github.com/opensim-org/opensim-core)
2. The opensim-reflex-controllers plugin (https://github.com/msdemers/opensim-reflex-controllers)
3. a git client
4. Source code (git clone https://github.com/msdemers/landing-optimization-framework)
5. [CMake](http://www.cmake.org/)

###Building Steps

####Linux (assuming you have gcc and/or clang)
1. build/install OpenSim. I'll call the install location ~/OpenSim_Install
2. build/install the opensim-reflex-controllers plugin
3. Creat a directories for cloning and building the landing-optimization-framework sourcecode into binaries
```
$ mkdir ~/landingOptimization
$ cd ~/landingOptimization
```
4. If you haven't yet, clone this repository, for example to ~/landingOptimization/source
```
$ git clone https://github.com/msdemers/landing-optimization-framework ~/landingOptimization/source
```
you can optionally check out a specific branch
```
$ git checkout [branch name]
```
5. Create a place to build the optimization binaries
```
$ mkdir build
```
6. Run CMake, choosing the IDE project system of your choice, or simply UNIX Make Files will serve fine.
  1. point at ~/landingOptimization/source for the source directory. This is where the high-level CMakeLists.txt lives
  2. point at ~/landingOptimization/build for the build directory
  3. Check the CMAKE_INSTALL_PREFIX variable.  This is where your the Install command will copy the binaries after they have been built.
  4. push *Configure* untill nothing appears red
  5. push *Generate*
  
7. Navigate to the build directory
```
$ cd build
```
8. Build the project. This may mean opening your IDE or using your favorite command line build tools.  If using make files.
```
make install
```

