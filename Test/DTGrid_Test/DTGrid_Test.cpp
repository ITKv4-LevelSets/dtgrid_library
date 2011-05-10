/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/
#include <string>

#include "DTGrid.h"
#include "LevelSet3D.h"
#include "MeanCurvatureFlowScalarField3D.h"
#include "TrilinearInterpolator.h"




/** The same as Grids::DTGridTraitsDefault except that this Trait allows processing without a safe band (includeSafeBand=false).
  * Note that this slows down stencil iteration somewhat and should not be enabled if not needed.
  * Processing without a safe band is required if rebuild() has not been called on the grid before iterating over it using stencil iterators.
  */
template<typename Data, typename Index, typename Real, typename UInt>
struct DTGridTraitsClosedLevelSet
{
    typedef Data DataType;
    typedef Index IndexType;
    typedef Real RealType;
    typedef Data *DataPtr;
    typedef UInt UIntType;

    /** The width of dilation used when rebuilding the narrow band (corresponds to the max number of voxels the interface can move between rebuilds) */
    static const UInt rebuildDilationWidth = 1;
    /** If the region includes an extra band of voxels around the gamma tube (added by a narrowBandRebuild() operation), some optimizations can be applied for iterators that rely on the DTGrid storing a signed distance field, and includeSafeBand should be defined. Note that when using the rebuild method supplied with the DT-Grid, the narrow band will always include a safeband! Furthermore note that when using FMM-based reinitialization, includeSafeBand must be disabled because FMM will not add a safeband. */
    static const bool includeSafeBand = false;
    /** If this is defined, the calls to increment operations in the updateStencil() method are inlined directly in the updateStencilMethod. Improves performance slightly. */
    static const bool inlineUpdateStencilCalls = true;
    /** The maximum number of grid points allowed in a dilation operation. Used to allocate static arrays. */
    static const int maxDilationWidth = 10;
    static const Grids::DTGridSearchType randomAccessType = Grids::USE_LINEAR_SEARCH;
    //static const Grids::DTGridSearchType randomAccessType = Grids::USE_BINARY_SEARCH;
    /** If this is false, linear array indices in off-center stencil iterators are not maintained (and hence cannot be used). Speeds up iteration */
    static const bool maintainArrayIndexInStencil = false;
    /** The width in grid cells of the zero crossing band */
    static const int zeroCrossingWidth = 1;
    /** If this boolean is true, open level sets are supported in random access.
     *  Open level sets are always supported for iteration with a stencil iterator. */
    static const bool openLevelSetsSupport = false;
};



/** The same as Grids::DTGridTraitsDefault except that this Trait allows open level sets (openLevelSetsSupport=true)
  * and processing without a safe band (includeSafeBand=false). Note that both of these properties slow down processing somewhat.
  * In particular they should not be enabled if not needed.
  */
template<typename Data, typename Index, typename Real, typename UInt>
struct DTGridTraitsOpenLevelSet
{
    typedef Data DataType;
    typedef Index IndexType;
    typedef Real RealType;
    typedef Data *DataPtr;
    typedef UInt UIntType;

    /** The width of dilation used when rebuilding the narrow band (corresponds to the max number of voxels the interface can move between rebuilds) */
    static const UInt rebuildDilationWidth = 1;
    /** If the region includes an extra band of voxels around the gamma tube (added by a narrowBandRebuild() operation), some optimizations can be applied for iterators that rely on the DTGrid storing a signed distance field, and includeSafeBand should be defined. Note that when using the rebuild method supplied with the DT-Grid, the narrow band will always include a safeband! Furthermore note that when using FMM-based reinitialization, includeSafeBand must be disabled because FMM will not add a safeband. */
    static const bool includeSafeBand = false;
    /** If this is defined, the calls to increment operations in the updateStencil() method are inlined directly in the updateStencilMethod. Improves performance slightly. */
    static const bool inlineUpdateStencilCalls = true;
    /** The maximum number of grid points allowed in a dilation operation. Used to allocate static arrays. */
    static const int maxDilationWidth = 10;
    static const Grids::DTGridSearchType randomAccessType = Grids::USE_LINEAR_SEARCH;
    //static const Grids::DTGridSearchType randomAccessType = Grids::USE_BINARY_SEARCH;
    /** If this is false, linear array indices in off-center stencil iterators are not maintained (and hence cannot be used). Speeds up iteration */
    static const bool maintainArrayIndexInStencil = false;
    /** The width in grid cells of the zero crossing band */
    static const int zeroCrossingWidth = 1;
    /** If this boolean is true, open level sets are supported in random access.
     *  Open level sets are always supported for iteration with a stencil iterator. */
    static const bool openLevelSetsSupport = true;
};

template< class TGrid >
void meanValue();

template< class TGrid >
void augmentedLevelSet();

template< class TGrid >
void randomAccess();

template< class TGrid >
void averageMeanCurvature();

template< class TGrid >
void localizedLevelSetComputation();

template< class TGrid >
void randomAccessOpenLevelSet();

template< class TGrid >
void rebuildOpenLevelSet();

template< class TGrid >
void csgOperation();

template< class TGrid >
void csgOperationWithTransform();

template< class TGrid >
void csgOperationOpenLevelSet();

template< class TGrid >
void neighborAccessOpenLevelSet();

template< class TGrid >
void closestPoint();


static string inputFileName;

/**
 *  Note that the methods that allow for open (unenclosed) level sets use the MyOpenGrid definition, whereas those methods that do not allow
 *  for open level sets uses the MyGrid definition. Open level sets are enabled through the DT-Grid's traits object
 *  (see DTGridTraitsClosedLevelSet and DTGridTraitsOpenLevelSet above), and open level sets are disabled per default as performance is better
 *  whenever it is disabled.
 */
int main(int argc, char* argv[])
{
  typedef float Real;
  typedef float Data;
  typedef short Index;
  typedef unsigned int UInt;

  ////////////////////////////////////////////////////////////////////////////////////////
  typedef DTGridTraitsOpenLevelSet<Data, Index, Real, UInt> MyTraitsOpenLevelSet;

  typedef Grids::DTGrid<MyTraitsOpenLevelSet> MyOpenGrid;
//  typedef MyOpenGrid::TubeIterator MyOpenTubeIterator;
//  typedef MyOpenGrid::EntireTubeIterator MyOpenEntireTubeIterator;
//  typedef MyOpenGrid::BetaTubeIterator MyOpenBetaTubeIterator;
//  typedef MyOpenGrid::GammaTubeIterator MyOpenGammaTubeIterator;
//  typedef MyOpenGrid::ZeroCrossingIterator MyOpenZeroCrossingIterator;
//  typedef MyOpenGrid::InitParams MyOpenInitParams;

//  typedef LevelSet::LevelSet3D<MyOpenGrid> MyOpenLevelSet;
//  typedef MyOpenLevelSet::InitParams MyOpenLSInitParams;
//  typedef MyOpenLevelSet::SimulationParameters MyOpenSimulationParameters;
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  typedef DTGridTraitsClosedLevelSet<Data, Index, Real, UInt> MyTraitsClosedLevelSet;

  typedef Grids::DTGrid<MyTraitsClosedLevelSet> MyGrid;
//  typedef MyGrid::TubeIterator MyTubeIterator;
//  typedef MyGrid::EntireTubeIterator MyEntireTubeIterator;
//  typedef MyGrid::BetaTubeIterator MyBetaTubeIterator;
//  typedef MyGrid::GammaTubeIterator MyGammaTubeIterator;
//  typedef MyGrid::ZeroCrossingIterator MyZeroCrossingIterator;
//  typedef MyGrid::InitParams MyInitParams;

//  typedef LevelSet::LevelSet3D<MyGrid> MyLevelSet;
//  typedef MyLevelSet::InitParams MyLSInitParams;
//  typedef MyLevelSet::SimulationParameters MySimulationParameters;
  ////////////////////////////////////////////////////////////////////////////////////////

    int testType;

    if (argc != 3)
    {
  std::cout << "Usage: DTGrid_Test <svol-file-name> <test-id>" << std::endl;
  std::cout << "Tests: " << std::endl;
  std::cout << " 1: Iteration (compute mean value)" << std::endl;
  std::cout << " 2: Stencil Iteration (compute average mean curvature)" << std::endl;
  std::cout << " 3: Random Access (compute mean value)" << std::endl;
  std::cout << " 4: Augmented Level Set Computations" << std::endl;
  std::cout << " 5: Localized Level Set Computations on Open Level Sets" << std::endl;
  std::cout << " 6: Random Access in Open Level Sets" << std::endl;
  std::cout << " 7: Narrow Band Rebuild in Open Level Sets" << std::endl;
  std::cout << " 8: CSG Operations" << std::endl;
  std::cout << " 9: CSG Operations that take Transformations into Account" << std::endl;
  std::cout << "10: CSG Operations on Open Level Sets" << std::endl;
  std::cout << "11: Neighbor Access in Open Level Sets" << std::endl;
  std::cout << "12: Closest Point" << std::endl;
  std::cout << "13: Perform Unit Tests" << std::endl;
  exit(-1);
    }

    inputFileName = argv[1];
    testType = atoi(argv[2]);

    switch (testType)
    {
    case 1:
  meanValue<MyOpenGrid>();
  break;
    case 2:
  averageMeanCurvature<MyOpenGrid>();
  break;
    case 3:
  randomAccess<MyOpenGrid>();
  break;
    case 4:
  augmentedLevelSet<MyOpenGrid>();
  break;
    case 5:
  localizedLevelSetComputation<MyOpenGrid>();
  break;
    case 6:
  randomAccessOpenLevelSet<MyOpenGrid>();
  break;
    case 7:
  rebuildOpenLevelSet<MyOpenGrid>();
  break;
    case 8:
  csgOperation<MyGrid>();
  break;
    case 9:
  csgOperationWithTransform<MyGrid>();
  break;
    case 10:
  csgOperationOpenLevelSet<MyOpenGrid>();
  break;
    case 11:
  neighborAccessOpenLevelSet<MyOpenGrid>();
  break;
    case 12:
  closestPoint<MyGrid>();
  break;
    case 13:
  meanValue<MyOpenGrid>();
  averageMeanCurvature<MyOpenGrid>();
  randomAccess<MyOpenGrid>();
  augmentedLevelSet<MyOpenGrid>();
  localizedLevelSetComputation<MyOpenGrid>();
  randomAccessOpenLevelSet<MyOpenGrid>();
  rebuildOpenLevelSet<MyOpenGrid>();
  csgOperation<MyGrid>();
  csgOperationWithTransform<MyGrid>();
  csgOperationOpenLevelSet<MyOpenGrid>();
  neighborAccessOpenLevelSet<MyOpenGrid>();
  closestPoint<MyGrid>();
  break;
    }

    return 1;
}



template< class TGrid >
void rebuildOpenLevelSet()
{
    typename TGrid::InitParams initParams;
    typename TGrid::Index bbox[3][2];
    typename TGrid::Index bboxDim[3];
    typename TGrid::UInt numValuesInOriginalSubset;
    typename TGrid::UInt numValuesInRebuiltSubset;


    TGrid *grid = new TGrid(initParams);
    grid->loadSparseVolume(inputFileName);

    // create open level set subset narrow band
    grid->boundingBox(bbox);
    grid->boundingBoxDim(bboxDim);
    // CASE 1
    //bbox[1][0] = bbox[1][0] + 10;
    //bbox[1][1] = bbox[1][0] + bboxDim[1]/4;
    //bbox[2][0] = bbox[2][0] + 2*bboxDim[2]/4;
    //bbox[2][1] = bbox[2][1] - bboxDim[2]/4;
    // CASE 2
    bbox[1][1] = bbox[1][0] + bboxDim[1]/4;
    // CASE 3
    //bbox[2][1] = bbox[2][0] + (bbox[2][1]-bbox[2][0])/4;
    TGrid* subset = grid->copy(bbox);

    subset = grid->copy(bbox);
    // set the open level set bbox on the subset.
    // this ensures that the level set does not expand beyond
    // this bounding box during dilation and rebuild operations
    subset->setOpenLevelSetBBox(bbox);

    // save subset as sparse volume for inspection
    subset->saveSparseVolume("rebuildopenlevelset_subset.svol");

    numValuesInOriginalSubset = subset->getNumValues();
    std::cout << "Number of values in original subset = " << numValuesInOriginalSubset << std::endl;
    // The assumption of the rebuildNarrowBandMethod is that the narrow band was fully contained
    // within the openLevelSetBoundingBox prior to the rebuild
    subset->rebuildNarrowBand();
    numValuesInRebuiltSubset = subset->getNumValues();
    std::cout << "Number of values in rebuilt subset = "
              << numValuesInRebuiltSubset
              << std::endl;

    subset->saveSparseVolume("rebuildopenlevelset_subset_rebuilt.svol");

    if (numValuesInRebuiltSubset < numValuesInOriginalSubset)
      {
      std::cout << "ERROR in rebuildOpenLevelSet():" << std::endl;
      std::cout << " Number of values in original subset <  Number of values in rebuilt subset" << std::endl;
      }
    else
      {
      std::cout << "No runtime errors in rebuildOpenLevelSet(), check result!"
                << std::endl;
      }

    delete grid;
    delete subset;
}



template< class TGrid >
void randomAccessOpenLevelSet()
{
    TGrid *grid, *subset;
    typename TGrid::InitParams initParams;
    typename TGrid::Index bbox[3][2];
    typename TGrid::Index bboxDim[3];
    typename TGrid::Index x, y, z;
    bool error = false;

    grid = new TGrid(initParams);
    grid->loadSparseVolume(inputFileName);

    // create open level set subset narrow band
    grid->boundingBox(bbox);
    grid->boundingBoxDim(bboxDim);
    // CASE 1
    bbox[1][0] = bbox[1][0] + 10;
    bbox[1][1] = bbox[1][0] + bboxDim[1]/4;
    bbox[2][0] = bbox[2][0] + 2*bboxDim[2]/4;
    bbox[2][1] = bbox[2][1] - bboxDim[2]/4;
    // CASE 2
    //bbox[1][1] = bbox[1][0] + bboxDim[1]/4;

    subset = grid->copy(bbox);

    // save subset as sparse volume for inspection
    subset->saveSparseVolume("randomacccessopenlevelset_subset.svol");

    // Now traverse through all grid points in the bbox of the subset level set and check
    // that the correct signed value is returned
    for (x = bbox[0][0]; x<=bbox[0][1] && !error; x++)
    {
  for (y = bbox[1][0]; y<=bbox[1][1] && !error; y++)
  {
      for (z = bbox[2][0]; z<=bbox[2][1] && !error; z++)
      {
    typename TGrid::Real value1 = (*grid)(x,y,z);
    typename TGrid::Real value2 = (*subset)(x,y,z);  // random access
    typename TGrid::Real value3;
    (*subset)(x,y,z,&value3);  // random access, approach 2

		if (value1 != value2 || value1 != value3)
		{
				std::cout << "ERROR in randomAccessOpenLevelSet() at ("
									<< x << "," << y << "," << z << ")" << std::endl;
				std::cout << "value1 = " << value1 << std::endl;
				std::cout << "value2 = " << value2 << std::endl;
				std::cout << "value3 = " << value3 << std::endl;
				error = true;
		}
			}
	}
		}
		if (!error)
		{
		std::cout << "No errors in randomAccessOpenLevelSet()!" << std::endl;
		}

    delete grid;
    delete subset;
}








/**
 *  This example shows how to perform localized level set operations on the DT-Grid.
 *  Since the DT-Grid encoding is fully sequential a subset of the narrow band must
 *  first be extracted. Local computations are then performed on this subset. Within a
 *  transition region that lies within a specified distance of the simulation bounding box
 *  the solution blends between the original model and the modified subset to allow for the subset to
 *  be pasted into the original model without visual seams. The blending strategy is currently implemented
 *  in the storeAuxiliaryValue() method of the LevelSet class. Note that in the case of topological changes
 *  and large deformations on the modified subset, this simple approach does not work, but the basic
 *  functionality for more advanced approaches is present. Finally
 *  the modified subset is then pasted back into the original level set. A smooth blending
 *  region close to the boundary of the subset ensures that it blends in smoothly with
 *  the original level set.
 *  In this example we perform a smoothing on the lower quarter of a model.
 *  Note that since the resulting subset level set can be open (ie not form a closed surface)
 *  it is important that 'openLevelSetsSupport = true' in the DTGrid traits class.
 */
template< class TGrid >
void localizedLevelSetComputation()
{
    typedef LevelSet::MeanCurvatureFlowScalarField3D<
      typename TGrid::Real,
      typename TGrid::Index,
      TGrid,
      typename TGrid::GammaTubeIterator> MyScalarField;

    TGrid *grid, *subset;
    typename TGrid::InitParams initParams;
    typename TGrid::Index bbox[3][2];

    typedef LevelSet::LevelSet3D<TGrid> MyOpenLevelSet;
    typename MyOpenLevelSet::InitParams lsInitParams;
    MyOpenLevelSet *ls;
    typename MyOpenLevelSet::SimulationParameters simParameters;
    MyScalarField *sf;

    typedef typename TGrid::Real Real;


    grid = new TGrid(initParams);
    grid->loadSparseVolume(inputFileName);

    // create subset narrow band
    grid->boundingBox(bbox);
    bbox[2][1] = bbox[2][0] + (bbox[2][1]-bbox[2][0])/2;
    subset = grid->copy(bbox);
    // set the open level set bbox on the subset.
    // this ensures that the level set does not expand beyond
    // this bounding box during dilation and rebuild operations
    subset->setOpenLevelSetBBox(bbox);

    // perform level set computations on the narrow band subset
    lsInitParams.numReinitIterations = 3;
    lsInitParams.useVanillaPeng = false;
    lsInitParams.reinitCFLNumber = static_cast<Real>(0.3);
    lsInitParams.propagateCFLNumber = static_cast<Real>(0.9);
    lsInitParams.verbose = true;
    lsInitParams.maxAbsDeterminationMethod = MyOpenLevelSet::MAXABS_NORM2;

    simParameters.numIterations = 2;
    simParameters.saveIntermediateSvols = true;
    simParameters.dt = static_cast<Real>(0.01);
    simParameters.maxTime = std::numeric_limits<Real>::max();
    simParameters.useSimBBox = true;
    simParameters.bboxSmoothTransitionWidth[0] = 0;
    simParameters.bboxSmoothTransitionWidth[1] = 10;   // smooth transition of 10 grid cells. I choose this rather conservatively
    simParameters.outputBaseFileName = "localizedLevelSetComputation";
    // We only want the smooth transition to be significant along the
    // border where we cut the model, therefore we expand the simBBox along the other directions
    // to avoid that features close to the other simBBox boundaries will not be smoothed.
    // Note that this is a very simple strategy and ideally a mask identifying the cut exactly should
    // be used to determine the transition region.
    simParameters.simBBox[0][0] = bbox[0][0]-simParameters.bboxSmoothTransitionWidth[1];
    simParameters.simBBox[0][1] = bbox[0][1]+simParameters.bboxSmoothTransitionWidth[1];
    simParameters.simBBox[1][0] = bbox[1][0]-simParameters.bboxSmoothTransitionWidth[1];
    simParameters.simBBox[1][1] = bbox[1][1]+simParameters.bboxSmoothTransitionWidth[1];
    simParameters.simBBox[2][0] = bbox[2][0]-simParameters.bboxSmoothTransitionWidth[1];
    simParameters.simBBox[2][1] = bbox[2][1];

    sf = new MyScalarField();
    sf->init(subset);

    // No upwind differences are used here, always second order central difference
    // SD_FIRSTORDER must be specified anyway, and indicates the stencil used for upwind differences, which as said above are not used here.
    lsInitParams.tdFormat = MyOpenLevelSet::TD_EULER;
    lsInitParams.sdFormat = MyOpenLevelSet::SD_FIRSTORDER;
    ls = new MyOpenLevelSet(subset, lsInitParams, false);  // false tells MyLevelSet that it does not own (and hence should not destruct) the subset grid
    ls->template propagateDt<
        Grids::SF_FIRSTORDER_CURVATURE,
        Grids::SF_FIRSTORDER,
        MyScalarField >(sf, simParameters);

    delete ls;
    delete sf;

    // paste the subset narrow band back into the original level set
    grid->paste(subset, bbox);
    delete subset;

    // save the final result as a sparse volume (svol) file
    grid->saveSparseVolume("localizedLevelSetComputation_final_result.svol");
    delete grid;

    std::cout << "No runtime errors in localizedLevelSetComputation(), check result!" << std::endl;
}






/**
 * This example demonstrates iterating over the narrow band to compute the mean value
 */
template< class TGrid >
void meanValue()
{
    typename TGrid::InitParams initParams;
    typename TGrid::Real mean = 0;

    TGrid* grid = new TGrid(initParams);
    grid->loadSparseVolume(inputFileName);

    typename TGrid::TubeIterator iter = grid->beginTubeIterator();
    typename TGrid::TubeIterator iend = grid->endTubeIterator();

    while (iter != iend)
    {
    typename TGrid::Real value = iter.getValue();
    mean += value;
    ++iter;
    }

    mean /= grid->getNumValues();

    std::cout << "mean = " << mean << std::endl;

    std::cout << "No runtime errors in meanValue(), check result!" << std::endl;
}


/**
 * This example demonstrates the random access operator to compute the mean value
 * and compares the result to the mean value computed using iteration
 */
template< class TGrid >
void randomAccess()
{
    TGrid *grid;
    typename TGrid::InitParams initParams;
    typename TGrid::TubeIterator iter, iend;
    typename TGrid::Real value, value2;
    typename TGrid::Real meanIter = 0;
    typename TGrid::Real meanRandom = 0;
    typename TGrid::Index x=0, y=0, z=0;
    bool error = false;

    grid = new TGrid(initParams);
    grid->loadSparseVolume(inputFileName);

    iter = grid->beginTubeIterator();
    iend = grid->endTubeIterator();

    while (iter != iend && !error)
    {
  // mean value by iteration
  value = iter.getValue();
  meanIter += value;

	// mean value by random access
	iter.getIndex(&x,&y,&z);
	value = (*grid)(x, y, z);  // random access
	(*grid)(x, y, z, &value2); // random access, approach 2
	if (value != value2)
	{
			error = true;
	}

	meanRandom += value;

  iter++;
    }
    if (error)
    {
  std::cout << "ERROR in randomAccess(), value = " << value << ", value2 = " << value2 << " at (" << x << "," << y << "," << z << ")" << std::endl;
    }
    else
    {
  if (meanIter != meanRandom)
  {
      std::cout << "ERROR in randomAccess(), meanIter != meanRandom" << std::endl;
      std::cout << "meanIter = " << meanIter << std::endl;
      std::cout << "meanRandom = " << meanRandom << std::endl;
  }
  else
  {
      std::cout << "No error in randomAccess()!" << std::endl;
  }
    }

    meanIter /= grid->getNumValues();

    std::cout << "meanIter = " << meanIter << std::endl;
}



/**
 *  Augmenting a level set with additional data is useful in many applications.
 *  Augmented data can be e.g. texture, colors, normals or velocities from a fluid simulation.
 *  This example shows how the DT-Grid can be augmented with additional data arrays and how to maintain
 *  this data during narrow band rebuild.
 *  In this simple example, the data is simple an array of 0 on the surface, +1 outside the interface and -1 inside the interface. */
template< class TGrid >
void augmentedLevelSet()
{
    TGrid *grid;
    typename TGrid::InitParams initParams;
    typename TGrid::UInt size, sizeAfterRebuild;
    typename TGrid::UInt *permutationArray;
    int *augmentedData, *augmentedDataAfterRebuild;
    bool error = false;

    grid = new TGrid(initParams);
    grid->loadSparseVolume(inputFileName);

    size = grid->getNumValues();
    augmentedData = new int[size+2];  // make room for inside/outside values
    for (typename TGrid::UInt i=0; i<size; i++)
    {
  augmentedData[i] = ( grid->getValue(i) < 0 ? -1 : ( grid->getValue(i) > 0 ? 1 : 0) );
    }
    // set inside value
    augmentedData[size] = -1;
    // set outside value
    augmentedData[size+1] = 1;


    // After the rebuild, the permutation array will contain
    // an entry for each grid point in the new narrow band.
    // Each entry contains an index into the old augmented data array.
    // If a grid point in the new narrow band did not exist in the old
    // narrow band its value is 'size' for an inside value, and 'size+1'
    // an outside value.
    grid->rebuildNarrowBand(&permutationArray);


    // construct augmented data after rebuild
    sizeAfterRebuild = grid->getNumValues();
    augmentedDataAfterRebuild = new int[sizeAfterRebuild + 2];
    for (typename TGrid::UInt i=0; i<sizeAfterRebuild; i++)
    {
  augmentedDataAfterRebuild[i] = augmentedData[ permutationArray[i] ];
    }
    // set inside value
    augmentedDataAfterRebuild[sizeAfterRebuild] = -1;
    // set outside value
    augmentedDataAfterRebuild[sizeAfterRebuild+1] = 1;

    // Test
    for (typename TGrid::UInt i=0; i<sizeAfterRebuild && !error; i++)
    {
  if (permutationArray[i]>size+2)
  {
      std::cout << "ERROR in augmentedLevelSet() test 1 at i = " << i << std::endl;
      std::cout << "permutation array value = " << permutationArray[i] << std::endl;
      error = true;
  }
    }
    // Test
    for (typename TGrid::UInt i=0; i<sizeAfterRebuild && !error; i++)
    {
  if (grid->getValue(i) < 0 && augmentedDataAfterRebuild[i]>=0 ||
      grid->getValue(i) > 0 && augmentedDataAfterRebuild[i]<=0 ||
      grid->getValue(i) == 0 && augmentedDataAfterRebuild[i]!=0)
  {
      std::cout << "ERROR in augmentedLevelSet() test 2 at i = " << i << std::endl;
      std::cout << "Grid value = " << grid->getValue(i) << std::endl;
      std::cout << "Augmented data value = " << augmentedDataAfterRebuild[i] << std::endl;
      std::cout << "sizeAfterRebuild = " << sizeAfterRebuild << std::endl;
      std::cout << "permutation array value = " << permutationArray[i] << std::endl;
      error = true;
  }
    }
    if (!error)
    {
  std::cout << "No errors in augmentedLevelSet()!" << std::endl;
    }


    delete[] augmentedData;
    delete[] permutationArray;
    delete grid;
    delete[] augmentedDataAfterRebuild;
}






/** This example demonstrates iteration with a stencil over the narrow band
 *  to compute the average mean curvature */
template< class TGrid >
void averageMeanCurvature()
{
    TGrid *grid;
    typename TGrid::EntireTubeIterator iter, iend;
    typename TGrid::InitParams initParams;

    grid = new TGrid(initParams);
    grid->loadSparseVolume(inputFileName);

    // compute average mean curvature
    typename TGrid::Real dx = grid->getScale();
    typename TGrid::Real epsilon = (typename TGrid::Real)(1.5 * dx);
    typename TGrid::Real volElement = dx * dx * dx;
    typename TGrid::Real deltaVal;
    typename TGrid::Real c1 = (typename TGrid::Real)(1.0 /
                              (typename TGrid::Real(2.0) * epsilon));
    typename TGrid::Real c2 = (typename TGrid::Real)(M_PI / epsilon);
    typename TGrid::Real H;
    typename TGrid::Real gradLen;
    typename TGrid::Real area;
    typename TGrid::Real cellArea;
    typename TGrid::Real phiVal;
    typename TGrid::Real averageMeanCurvature;

    iter = grid->template beginEntireTube<Grids::SF_FIRSTORDER_CURVATURE, false>();
    iend = grid->template endEntireTube<Grids::SF_FIRSTORDER_CURVATURE, false>();

    averageMeanCurvature = 0;
    area = 0;

    while (iter != iend)
    {
  phiVal = iter.getValue();

	if (fabs(phiVal) <= epsilon)
	{
			deltaVal = c1 + c1 * cos( phiVal * c2 );
			H = iter.template meanCurvature<Grids::SF_FIRSTORDER_CURVATURE>();
			gradLen = iter.gradientLength();
			cellArea = deltaVal * gradLen * volElement;
			area += cellArea;
			averageMeanCurvature +=  H * cellArea;
	}

  iter.template operator++<Grids::SF_FIRSTORDER_CURVATURE, false>();
    }

    averageMeanCurvature /= area;

    std::cout << "Average mean Curvature is " << averageMeanCurvature << std::endl;

    std::cout << "No runtime errors in averageMeanCurvature(), check result!" << std::endl;

    delete grid;
}


/** For the purpose of demonstrating CSG operations a level set is loaded
 *  into two separate DT-Grids. The second DT-Grid is translated by changing the
 *  coordinates of its grid points and used
 *  to perform CSG operations union, difference and intersection with the
 *  first DT-Grid.
 *  The CSG operation performed here (CSG_GA, where GA stands for Grid Aligned)
 *  does not take transformations into account and performs the CSG operations
 *  on the basis of the voxel coordinates and their values alone. This method is
 *  fast because the grids live in the same coordinate system. The lexicographic
 *  storage order can hen be exploited and the CSG operation
 *  can be performed in a time complexity that is linear in the number of grid points
 *  in the two DT-Grids (hence the time complexity scales with the interface size).
 *  Note that this method does not recompute the signed distance function!
 */
template< class TGrid >
void csgOperation()
{
    TGrid *grid1, *grid2;
    typename TGrid::EntireTubeIterator iter, iend;
    typename TGrid::InitParams initParams;
    typename TGrid::Index bboxDim[3];

    grid1 = new TGrid(initParams);
    grid1->loadSparseVolume(inputFileName);
    grid1->boundingBoxDim(bboxDim);
    bboxDim[0] /= 8;
    bboxDim[1] /= 8;
    bboxDim[2] /= 8;
    grid1->translateGrid(bboxDim);

    // Union
    grid2 = new TGrid(initParams);
    grid2->loadSparseVolume(inputFileName);
    std::cout << "Union.."; std::cout.flush();
    grid2->CSG_GA(*grid1, TGrid::CSG_UNION);
    std::cout << "done" << std::endl;
    grid2->saveSparseVolume("union_result.svol");
    delete grid2;

    // Intersection
    grid2 = new TGrid(initParams);
    grid2->loadSparseVolume(inputFileName);
    std::cout << "Intersection.."; std::cout.flush();
    grid2->CSG_GA(*grid1, TGrid::CSG_INTERSECTION);
    std::cout << "done" << std::endl;
    grid2->saveSparseVolume("intersection_result.svol");
    delete grid2;

    // Difference
    grid2 = new TGrid(initParams);
    grid2->loadSparseVolume(inputFileName);
    std::cout << "Difference.."; std::cout.flush();
    grid2->CSG_GA(*grid1, TGrid::CSG_DIFFERENCE);
    std::cout << "done" << std::endl;
    grid2->saveSparseVolume("difference_result.svol");
    delete grid2;

    delete grid1;

    std::cout << "No runtime errors in csgOperation(), check result!" << std::endl;
}




/** For the purpose of demonstrating CSG operations involving transformations,
 *  a level set is loaded into two separate DT-Grids. The second DT-Grid is rotated and used
 *  to perform CSG operations union, difference and intersection with the
 *  first DT-Grid.
 *  The CSG operation performed here is not as fast as CSG_GA because it does not
 *  assume that the two DT-Grids live in the same coordinate system. Hence the lexicographic order
 *  can only be exploited in one of the DT-Grids. For the other grid we have to traverse the entire bounding box.
 *  Hence when doing these type of CSG operations with transformations it is an advantage to use the grid with
 *  the smallest bounding box as the argument to the CSG method.
 *  Note that this method does not recompute the signed distance function!
 */
template< class TGrid >
void csgOperationWithTransform()
{
    TGrid *grid1, *grid2;
    typename TGrid::EntireTubeIterator iter, iend;
    typename TGrid::InitParams initParams;

    typedef typename TGrid::Real Real;
    typename Math::TrilinearInterpolator<Real, Real> triInterp;

    grid1 = new TGrid(initParams);
    grid1->loadSparseVolume(inputFileName);
    grid1->setRotationXYZ(1,1,1);  // angles specified in radians

    // Union
    grid2 = new TGrid(initParams);
    grid2->loadSparseVolume(inputFileName);
    std::cout << "Union.."; std::cout.flush();
    grid2->CSG(*grid1, TGrid::CSG_UNION, &triInterp);
    std::cout << "done" << std::endl;
    grid2->saveSparseVolume("union_transform_result.svol");
    delete grid2;

    // Intersection
    grid2 = new TGrid(initParams);
    grid2->loadSparseVolume(inputFileName);
    std::cout << "Intersection.."; std::cout.flush();
    grid2->CSG(*grid1, TGrid::CSG_INTERSECTION, &triInterp);
    std::cout << "done" << std::endl;
    grid2->saveSparseVolume("intersection_transform_result.svol");
    delete grid2;

    // Difference
    grid2 = new TGrid(initParams);
    grid2->loadSparseVolume(inputFileName);
    std::cout << "Difference.."; std::cout.flush();
    grid2->CSG(*grid1, TGrid::CSG_DIFFERENCE, &triInterp);
    std::cout << "done" << std::endl;
    grid2->saveSparseVolume("difference_transform_result.svol");
    delete grid2;

    delete grid1;

    std::cout << "No runtime errors in csgOperationWithTransform(), check result!" << std::endl;
}





/** For the purpose of demonstrating CSG operations on an open level set, a level set is loaded
 *  into two separate DT-Grids. The second DT-Grid is translated by changing the
 *  coordinates of its grid points and used
 *  to perform CSG operations union, difference and intersection with the
 *  first DT-Grid.
 *  The CSG operation performed here (CSG_GA, where GA stands for Grid Aligned)
 *  does not take transformations into account and performs the CSG operations
 *  on the basis of the voxel coordinates and their values alone. This method is
 *  relatively fast because the grids live in the same coordinate system. However it is not as fast
 *  as CSG operations using closed level sets, because several different cases need to be taken
 *  into account. The lexicographic storage order can hen be exploited and the CSG operation
 *  can be performed in a time complexity that is linear in the number of grid points
 *  in the two DT-Grids (hence the time complexity scales with the interface size).
 *  Note that this method does not recompute the signed distance function!
 */
template< class TGrid >
void csgOperationOpenLevelSet()
{
    TGrid *grid1, *grid2, *subset1, *subset2, *subset3;
    typename TGrid::TubeIterator iter, iend, iter2, iend2;
    typename TGrid::InitParams initParams;
    typename TGrid::Index bboxDim[3];
    typename TGrid::Index bbox[3][2];

    grid1 = new TGrid(initParams);
    grid1->loadSparseVolume(inputFileName);
    grid1->boundingBoxDim(bboxDim);
    bboxDim[0] /= 8;
    bboxDim[1] /= 8;
    bboxDim[2] /= 8;
    grid1->translateGrid(bboxDim);

    // create subset narrow band (open level set) of grid1
    // and set the open level set bbox on the subset.
    grid1->boundingBox(bbox);
    bbox[2][1] = bbox[2][0] + (bbox[2][1]-bbox[2][0])/2;
    bbox[2][0] = 0;
    subset1 = grid1->copy(bbox);
    subset1->setOpenLevelSetBBox(bbox);
    subset1->saveSparseVolume("subset1.svol");

    // Union
    grid2 = new TGrid(initParams);
    grid2->loadSparseVolume(inputFileName);
    subset2 = grid2->copy(bbox);
    subset2->setOpenLevelSetBBox(bbox);
    std::cout << "Union.."; std::cout.flush();
    subset2->saveSparseVolume("subset2.svol");
    subset2->CSG_GA(*subset1, TGrid::CSG_UNION);
    std::cout << "done" << std::endl;
    subset2->saveSparseVolume("union_openls_result.svol");

    // Intersection
    grid2 = new TGrid(initParams);
    grid2->loadSparseVolume(inputFileName);
    subset2 = grid2->copy(bbox);
    subset2->setOpenLevelSetBBox(bbox);
    std::cout << "Intersection.."; std::cout.flush();
    subset2->CSG_GA(*subset1, TGrid::CSG_INTERSECTION);
    std::cout << "done" << std::endl;
    subset2->saveSparseVolume("intersection_openls_result.svol");
    delete grid2;
    delete subset2;

    // Difference
    grid2 = new TGrid(initParams);
    grid2->loadSparseVolume(inputFileName);
    subset2 = grid2->copy(bbox);
    subset2->setOpenLevelSetBBox(bbox);
    std::cout << "Difference.."; std::cout.flush();
    subset2->CSG_GA(*subset1, TGrid::CSG_DIFFERENCE);
    std::cout << "done" << std::endl;
    subset2->saveSparseVolume("difference_openls_result.svol");


    // Compare to difference with closed level set followed by a copy operation.
    // In particular this should give the same result.
    std::cout << "Difference with closed level set.."; std::cout.flush();
    grid2->CSG_GA(*grid1, TGrid::CSG_DIFFERENCE);
    std::cout << "done" << std::endl;
    subset3 = grid2->copy(bbox);

    iter = subset2->beginTubeIterator();
    iend = subset2->endTubeIterator();
    iter2 = subset3->beginTubeIterator();
    iend2 = subset3->endTubeIterator();

    int numPointsOk = 0;
    bool error = false;

    while ((iter!=iend || iter2!=iend2) && !error)
    {
  typename TGrid::Index c[3], c2[3];

	iter.getIndex(&c[0], &c[1], &c[2]);
	iter2.getIndex(&c2[0], &c2[1], &c2[2]);

	if (c[0] != c2[0] || c[1] != c2[1] || c[2] != c2[2])
	{
			error = true;
			std::cout << "ERROR in csgOperationOpenLevelSet():" << std::endl;
			std::cout << "point from open level set: " << c[0] << " " << c[1] << " " << c[2] << std::endl;
			std::cout << "point from closed level set: " << c2[0] << " " << c2[1] << " " << c2[2] << std::endl;
			std::cout << "NumPointsOK = " << numPointsOk << std::endl;
			std::cout << std::endl;
	}
	else
	{
			++numPointsOk;
	}

	if (iter!=iend)
			++iter;
	if (iter2!=iend2)
			++iter2;
		}
		if (!error)
		{
	std::cout << "No errors in csgOperationOpenLevelSet()!" << std::endl;
		}

    delete grid2;
    delete subset2;
    delete subset3;

    delete grid1;
    delete subset1;
}



/** This example illustrates neighbor access into an open level set.
 *  Note that neighbor access is generally faster than random access, since information
 *  about the current grid point can be used to retrieve neighbor information. */
template< class TGrid >
void neighborAccessOpenLevelSet()
{
    TGrid *grid, *subset;
    typename TGrid::InitParams initParams;
    typename TGrid::Index bbox[3][2];
    typename TGrid::Index bboxDim[3];
    typename TGrid::Index x, y, z;
    typedef typename TGrid::Real Real;
    bool error = false;

    grid = new TGrid(initParams);
    grid->loadSparseVolume(inputFileName);

    // create open level set subset narrow band
    grid->boundingBox(bbox);
    grid->boundingBoxDim(bboxDim);
    // CASE 1
    bbox[1][0] = bbox[1][0] + 10;
    bbox[1][1] = bbox[1][0] + bboxDim[1]/4;
    bbox[2][0] = bbox[2][0] + 2*bboxDim[2]/4;
    bbox[2][1] = bbox[2][1] - bboxDim[2]/4;
    // CASE 2
    //bbox[1][1] = bbox[1][0] + bboxDim[1]/4;

    subset = grid->copy(bbox);

    // save subset as sparse volume for inspection
    subset->saveSparseVolume("neighboracccessopenlevelset_subset.svol");

    // Now traverse through all grid points in the bbox of the subset level set and check
    // that the correct signed value is returned
    // Note that for an open level set it only makes sense to look up values inside its bounding box
    for (x = bbox[0][0]+1; x<=(bbox[0][1]-1) && !error; x++)
    {
  for (y = bbox[1][0]+1; y<=(bbox[1][1]-1) && !error; y++)
  {
      for (z = bbox[2][0]+1; z<=(bbox[2][1]-1) && !error; z++)
      {
    // Is this grid point included in the narrow band?
    if (grid->doesElementExist(x,y,z))
    {
        // Yes, grid point is included in the narrow band
        typename TGrid::Locator loc;

				// Retrieve locator at this grid point using random access.
				// If an iterator is available, a locator can also be constructed from the iterator using the getLocator() method on the iterator.
				grid->getLocator(x, y, z, &loc);

				// get neighbor (x-1,y,z)
				Real xm = grid->getXM(x, y, z, loc.ic1D);
				// get neighbor (x+1,y,z)
				Real xp = grid->getXP(x, y, z, loc.ic1D);
				// get neighbor (x,y-1,z)
				Real ym = grid->getYM(y, z, loc.ic2D);
				// get neighbor (x,y+1,z)
				Real yp = grid->getYP(y, z, loc.ic2D);
				// get neighbor (x,y,z-1)
				Real zm = grid->getZM(z, loc.ic3D);
				// get neighbor (x,y,z+1)
				Real zp = grid->getZP(z, loc.ic3D);

				Real xm2 = (*grid)(x-1, y, z);
				Real xp2 = (*grid)(x+1, y, z);
				Real ym2 = (*grid)(x, y-1, z);
				Real yp2 = (*grid)(x, y+1, z);
				Real zm2 = (*grid)(x, y, z-1);
				Real zp2 = (*grid)(x, y, z+1);

				if (xm!=xm2)
				{
			std::cout << "ERROR in neighborAccessOpenLevelSet() at (" << x << "," << y << "," << z << ")" << std::endl;
			std::cout << "xm1 = " << xm << std::endl;
			std::cout << "xm2 = " << xm2 << std::endl;
			error = true;
				}
				if (xp!=xp2)
				{
			std::cout << "ERROR in neighborAccessOpenLevelSet() at (" << x << "," << y << "," << z << ")" << std::endl;
			std::cout << "xp1 = " << xp << std::endl;
			std::cout << "xp2 = " << xp2 << std::endl;
			error = true;
				}

				if (ym!=ym2)
				{
			std::cout << "ERROR in neighborAccessOpenLevelSet() at (" << x << "," << y << "," << z << ")" << std::endl;
			std::cout << "ym1 = " << ym << std::endl;
			std::cout << "ym2 = " << ym2 << std::endl;
			error = true;
				}
				if (yp!=yp2)
				{
			std::cout << "ERROR in neighborAccessOpenLevelSet() at (" << x << "," << y << "," << z << ")" << std::endl;
			std::cout << "yp1 = " << yp << std::endl;
			std::cout << "yp2 = " << yp2 << std::endl;
			error = true;
				}

				if (zm!=zm2)
				{
			std::cout << "ERROR in neighborAccessOpenLevelSet() at (" << x << "," << y << "," << z << ")" << std::endl;
			std::cout << "zm1 = " << zm << std::endl;
			std::cout << "zm2 = " << zm2 << std::endl;
			error = true;
				}
				if (zp!=zp2)
				{
			std::cout << "ERROR in neighborAccessOpenLevelSet() at (" << x << "," << y << "," << z << ")" << std::endl;
			std::cout << "zp1 = " << zp << std::endl;
			std::cout << "zp2 = " << zp2 << std::endl;
			error = true;
				}

		}
			}
	}
		}
		if (!error)
		{
	std::cout << "No errors in neighborAccessOpenLevelSet()!" << std::endl;
		}

    delete grid;
    delete subset;
}



template< class TGrid >
void closestPoint()
{
    TGrid *grid;
    typename TGrid::InitParams initParams;
    typename TGrid::Index bbox[3][2];
    typename TGrid::Index x, y, z;
    typename TGrid::Index p[3];               // grid point
    typename TGrid::Real n[3];               // normal vector at grid point
    typename TGrid::Real cd;                  // closest distance to interface at grid point
    bool inside;              // true if grid point is inside the object
    bool inNarrowBand;        // true if grid point is inside the narrow band

    typename TGrid::UInt numInNarrowBand = 0;
    typename TGrid::UInt numInsideObject = 0;

    grid = new TGrid(initParams);
    grid->loadSparseVolume(inputFileName);

    grid->boundingBox(bbox);

    for (x=bbox[0][0]; x<=bbox[0][1]; x++)
    {
  for (y=bbox[1][0]; y<=bbox[1][1]; y++)
  {
      for (z=bbox[2][0]; z<=bbox[2][1]; z++)
      {
    p[0] = x;
    p[1] = y;
    p[2] = z;
    grid->closestPoint(p, n, &cd, &inside, &inNarrowBand);

		if (inside)
				numInsideObject++;
		if (inNarrowBand)
				numInNarrowBand++;
			}
	}
		}

    std::cout << "numInsideObject = " << numInsideObject << std::endl;
    std::cout << "numInNarrowBand = " << numInNarrowBand << std::endl;

    if (numInNarrowBand != grid->getNumValues())
    {
  std::cout << "ERROR in closestPoint(), numInNarrowBand = " << numInNarrowBand << ", grid->getNumValues() = " << grid->getNumValues() << std::endl;
    }
    else
    {
  std::cout << "No errors in closestPoint()!" << std::endl;
    }

}
