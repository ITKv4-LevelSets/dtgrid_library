Introduction
============

This document is a brief guide to getting started using the DT-Grid library.
If you have any questions, comments, suggestions or encounter any bugs do not hesitate to email me (Michael Bang Nielsen at nielsenmb@gmail.com).

The DT-Grid is a data structure for high resolution level set simulations.
The memory usage is proportional to the size of the narrow band as opposed to the volume of the simulation domain,
and in addition the DT-Grid performs relatively fast due to cache coherency.
The DT-Grid and derivative work have been used in many applications ranging from surface processing to fluid simulation.
For more detail see for example the following papers and technical sketches:


Nielsen, M.B., Museth, K.: Dynamic Tubular Grid: An Efficient Data Structure and Algorithms for High Resolution Level Sets. Journal of Scientific Computing, Volume 26, Number 3, March 2006. Pages 261-299.

Houston, B., Nielsen, M.B., Batty, C., Nilsson, O., Museth, K.: Hierarchical RLE Level Set: A Compact and Versatile Deformable Surface Representation. ACM Transactions on Graphics, Volume 25, Number 1, January 2006. Pages 151-175.

Nielsen, M.B., Nilsson, O., Söderström, A., Museth, K.: Out-Of-Core and Compressed Level Set Simulations, ACM Transactions on Graphics 26(4), October 2007 (Presented at SIGGRAPH 2008, Partial Differential Equations paper session).

Nemitz, O., Nielsen, M. B., Rumpf, M., Whitaker, R.: Finite Element Methods on Very Large, Dynamic Tubular Grid Encoded Implicit Surfaces, SIAM Journal of Scientific Computing, Vol. 31, No. 3. pp. 2258-2281, 2009.

Nemitz, O., Nielsen, M. B., Rumpf, M., Whitaker, R.: Narrow Band Methods for PDEs on Very Large Implicit Surfaces. In Vision, Modeling and Visualization Proceedings, pages 171-180, November 2007.

Nielsen, M.B., Nilsson, O., Söderström, A., Museth, K.: Virtually Infinite Resolution Deformable Surfaces. In Proc. SIGGRAPH 2006 Sketches and Applications, August 2006.

Houston, B., Nielsen, M.B., Batty, C., Nilsson, O., Museth, K.: Gigantic Deformable Surfaces. In Proc. SIGGRAPH 2005 Sketches and Applications, August 2005.

Nielsen, M.B., Museth, K.: An Optimized, Grid Independent, Narrow Band Data Structure. In Proc. SIGRAD 2004 Sketches and Applications, November 2004.

Nielsen, M.B.: Efficient and High Resolution Level Set Simulations - Data Structures, Algorithms and Applications, PhD Dissertation, Aarhus University, December 2006.




Feature List
============

This package contains the following features:

- DT-Grid narrow band data structure
   - Support for iterators
   - Support for random and neighbor access
   - Support for narrow band dilation and rebuild
   - Support for first order one-sided, second order central and WENO finite difference schemes
   - Support for CSG operations
   - Support for open (ie unenclosed) level set volumes
     - Extration of open level set subvolumes from larger closed or open level set volumes
     - Random and neighbor access, csg and localized level set operations on open level set volumes
     - Insertion of open level set volumes back into the original level set volume
   - Support for augmented level sets
     - Additional fields such as texture coordinates, normals etc can be stored in additional to the level set function.
       These fields are automatically re-allocated whenever the narrow band changes.

- Tool for conversion of closed .obj and .ply mesh files to DT-Grids

- Tool for conversion of DT-Grids to .ply mesh files

- Tools that illustrate
  - Level set motion in the normal direction
    - Mean curvature flow
    - Volume conserving mean curvature flow
    - Shape metamorphosis
  - Level set advection in an externally generated velocity field
    - Enright Test

- Test program that illustrates the following features
  - Iterating over the narrow band
  - Iterating a stencil over the narrow band
  - Random and neighbor access (on open and closed level sets)
  - CSG operations (on open and closed level sets)
  - Computing the closest point on a surface
  - Level set computations (on open and closed level sets)
  - Localized level set computations on subsets of the level set narrow band (open level sets)
  - Rebuilding the narrow band (on augmented, open and closed level sets)




Upcoming Features
=================

I plan to add the following features in the future:

- More elaborate documentation
- Parallel implementation using Intel Thread Building Blocks
- Out-of-core scan converter that allows larger models to be scan converted
- Fast Marching implementation on the DT-Grid
- DT-Grid memory optimizations
- DT-Grid visualizer
- Compressed storage of DT-Grids
- Robust scan conversion of meshes with holes and self-intersections.




Compilation and Installation
============================

All downloads referred to below are available at the DT-Grid Google Code project.

The DT-Grid library has been successfully compiled and installed on the following platforms:
Linux: g++ (GCC) 4.1.2
Windows XP 32 bit and Windows Vista Ultimate 64 bit: Visual Studio 2008 version 9.0.
Mac OS X 10.5.8: g++ 4.0.1

In addition to the core DT-Grid library, to install, you will need to download the additional software packages:
CPT.zip
obj-0.1.zip
ply-0.1.zip
These software packages should be placed in the directory named External in the DT-Grid library source tree.

In addition you should download the source code for the DT-Grid based level set tools:
Tools.zip
MeshTools.zip
This source code should be placed in the root directory of the DT-Grid library source tree.

For getting started, you can also download a closed (watertight) version of the Stanford Bunny:
stanford_bunny.obj

After correctly unzipping these files, the root of your DT-Grid library source tree should look like this:

Algorithms/      Doxyfile   Test/                 development.txt
Core/            External/  Tools/                getting_started.txt
DataStructures/  Makefile   acknowledgements.txt  version_history.txt

It is also required to have the Boost library installed: www.boost.org
For Windows, a binary installer can be found at http://www.boostpro.com/download
Only header files are required by the DT-Grid library which simplifies the installation, see http://www.boost.org/doc/libs/1_40_0/more/getting_started/index.html


Compiling and Installing on Windows:
===================================

1) In case your Windows system does not include the C++ Technical Report 1 (TR1) Library Extensions (see http://en.wikipedia.org/wiki/C%2B%2B_Technical_Report_1), you can download the following feature pack from the official microsoft site that includes TR1:
http://www.microsoft.com/downloads/details.aspx?FamilyId=D466226B-8DAB-445F-A7B4-448B326C48E7&displaylang=en

2) In case you want the binaries to be copied to a common directory after compilation (post-build-event), define the environment variable DTGridBinariesPath to point to this directory. Visual Studio will then automatically copy binaries to this directory when compiling in Release mode. You can set environment variables by opening Control Panel -> System -> Advanced -> Environment Variables. Note that you will need to restart your Visual Studio in case it is already open.

3) Start Visual Studio and open the solution file Test/DTGrid_Test/DTGrid_Test.sln

4) Build the project. For debugging build in Debug mode, for performance build in Release mode.


Compiling and Installing on Linux:
=================================

1) In case you want the binaries to be copied to a common directory after compilation, define the environment variable DTGridBinariesPath to point to this directory. In a BASH shell this is done by typing:

DTGridBinariesPath=YOUR_DIRECTORY
export DTGridBinariesPath

2) Build the binaries. This is done by typing:

make all

in the root directory of the DT-Grid library source tree.

3) To copy the binaries to the directory specified by DTGridBinariesPath, type
make install


Compiling and Installing on Mac OS X:
=================================

1) Open a Terminal window.

2) In case you want the binaries to be copied to a common directory after compilation, define the environment variable DTGridBinariesPath to point to this directory. This is done by typing:

DTGridBinariesPath=YOUR_DIRECTORY
export DTGridBinariesPath

3) Unless the boost libraries are in your default include path, you will need to define the environment variable BOOST_ROOT. This is done by typing:

BOOST_ROOT=PATH_TO_BOOST_ROOT
export BOOST_ROOT

4) Build the binaries. This is done by typing:

make all

in the root directory of the DT-Grid library source tree.
NOTE: I did in some cases experience a problem when compiling the file Tools/MeshToDTGrid/PlyObjToDTGrid.cpp. The compiler exited with the error message "virtual memory exhausted: Cannot allocate memory". A remedy is to remove the -03 optimization flag in the makefile Tools/MeshToDTGrid/Makefile, although this will unfortunately degrade the performance of the mesh to DT-Grid converter.

5) To copy the binaries to the directory specified by DTGridBinariesPath, type:
make install




Building the Documentation
==========================

The DT-Grid library supports Doxygen (http://www.stack.nl/~dimitri/doxygen/) style comments and includes a Doxyfile.
To build the documentation for the library, type:

Doxygen Doxyfile

in the root directory of the DT-Grid library source tree.




Using the DT-Grid Library
========================================

To get started I suggest you proceed with the "Getting Started" section below. 
After that you should consult the source code of the DTGrid_Test program to see uses of iterators, random access, CSG operations and so on.
Finally the tools Erosion, EnrightTest etc show you how to setup a level set simulation.


Important note about the DT-Grid Traits object
==============================================

As is evident from the DTGrid_Test source code, the DT-Grid is templated by a traits object that you must define if not using the default. 
This traits object defines some of the basic features of the resulting class. It is important that you understand the basic properties of this traits object and you should consult the documentation in the DTGrid.h file or the generated Doxygen files. The default traits object (DTGridTraitsDefault defined in DTGrid.h) is defined for the best run-time performance. However you will have to define your own traits object in order to use e.g. open level sets or level sets that do not include a safe-band of grid points added by the rebuild method. See the DTGrid_Test program as well as the Erosion, EnrightTest etc. level set tools for how to setup and apply the traits object.




Getting Started
============

This brief introduction will guide you through a few steps to getting started with the DT-Grid library.
After having completed these steps you should be ready to begin using the library or developing your own applications, extensions and modifications.


1) Complete the Compilation and Installation procedure above and ensure that the executables are in your path.


2) Download a mesh that will be converted into a level set. E.g. the watertight Stanford Bunny supplied with this project.

   Alternatively I recommend you visit the Stanford Scanning Repository:
   http://graphics.stanford.edu/data/3Dscanrep/
   Note that the mesh you use must be closed with the current tools supplied in the DT-Grid library. If this is not the case, the MeshToDTGrid tool may cause an exception or generate corrupt data. Examples of models we have successfully converted are the Lucy, Asian Dragon, Thai Statue and Happy Buddha models.


3) Run the DT-Grid library unit tests to check if everything works as intended:
   First convert the mesh to a level set of some resolution (I use 150^3 below with a narrow band radius of 3):
   MeshToDTGrid.exe stanford_bunny.obj stanford_bunny_150 150 3

   Run the DT-Grid library unit tests:
   DTGrid_Test.exe stanford_bunny_150.svol 13

   To see which other tests are supported by the DTGrid_Test program, run it without any arguments. When you start using the library I recommend that you consult the source code of this test program to get confident in using iterators, random access etc.


4) Run the Enright Test:
   The Enright Test was proposed by Enright etal. in the paper "A Hybrid Particle Level Set Method for Improved Interface Capturing" (see http://physbam.stanford.edu/~fedkiw/papers/stanford2001-04.pdf).
   For an Enright Test at resolution 1024^3 see the video http://www.cs.au.dk/~bang/publications/enright_movie-200fps.mpeg

   Start the Enright test on a sphere in resolution 128^3 using the commandline:
   EnrightTest.exe sphere enright_result 128

   The EnrightTest.exe program produces a number of 'svol' files. An 'svol' file contains a DT-Grid.
   One way to visualize the DT-Grid is by converting it to a 'ply' mesh and subsequently loading it into a viewer.
   Alternatively it can be raytraced directly.

   To convert for example the file 'enright_result_iter_100_time_1.01.svol' to a 'ply' mesh, use the DTGridToMesh.exe application which is invoked as:
   DTGridToMesh.exe enright_result_iter_100_time_1.01.svol enright_result_iter_100_time_1.01
   This will generate a ply mesh file named enright_result_iter_100_time_1.01.ply which you can load into your preferred viewer.   

   To start the Enright Test on an object other than a sphere you need to first scan convert that object. For example you can convert the Stanford Bunny model to a mesh at resolution roughly 80^3 and narrow band radius 5 (here we use a radius of 5 as the Enright Test per default is run using a higher order numerical method) as follows:
   MeshToDTGrid.exe stanford_bunny.obj stanford_bunny_80 80 5

   To subsequently start the Enright Test at resolution 256^3 on this model issue the command:
   EnrightTest.exe stanford_bunny_80.svol bunny_enright_result 256

   An example of an older Enright Test we did at resolution 1024^3 with the Stanford Bunny can be found here:
   http://www.cs.au.dk/~bang/publications/BunnyEnright1024.avi


5) Run a shape metamorphosis:
   This example allows you to create a simple shape metamorphosis between two level sets. In the example provided here, the level sets used are both narrow band level sets meaning that the speed in-between the source and target shapes is constant. More interesting motion can be generated using the original technique proposed in the paper "A Level-Set Approach for the Metamorphosis of Solid Models" by Breen and Whitaker (http://www.cs.drexel.edu/~david/Papers/morph_tvcg.pdf).

   To morph between two differently sized Stanford Bunny models, first convert the Stanford Bunny mesh into two different level sets:
   MeshToDTGrid.exe stanford_bunny.obj stanford_bunny_100 100 3
   MeshToDTGrid.exe stanford_bunny.obj stanford_bunny_150 150 3

   To perform the shape metamorphosis, issue the command:
   Morph.exe stanford_bunny_100.svol stanford_bunny_150.svol stanford_bunny_morph 35


6) There are additional tools supplied in this distribution which you can check out:
   "Erosion" performs an erosion on a level set.
   "MeanCurvatureFlow" performs mean curvature flow (smoothing) on a level set.
   "VolumeConservingMeanCurvatureFlow" performs volume conserving mean curvature flow (volume conserving smoothing) on a level set.
