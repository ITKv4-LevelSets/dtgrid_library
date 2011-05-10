/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/
#include <string>

#include <DataStructures/Grids/DTGrid/HeaderFiles/DTGrid.h>
#include <Algorithms/LevelSet/HeaderFiles/LevelSet3D.h>
#include <Algorithms/LevelSet/ScalarFields/HeaderFiles/MeanCurvatureFlowScalarField3D.h>
#include <Algorithms/Math/Interpolators/HeaderFiles/TrilinearInterpolator.h>


using namespace std;

typedef float Real;
typedef float Data;
typedef short Index;
typedef unsigned int UInt;

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


////////////////////////////////////////////////////////////////////////////////////////
typedef DTGridTraitsOpenLevelSet<Data, Index, Real, UInt> MyTraitsOpenLevelSet;

typedef Grids::DTGrid<MyTraitsOpenLevelSet> MyOpenGrid;
typedef MyOpenGrid::TubeIterator MyOpenTubeIterator;
typedef MyOpenGrid::EntireTubeIterator MyOpenEntireTubeIterator;
typedef MyOpenGrid::BetaTubeIterator MyOpenBetaTubeIterator;
typedef MyOpenGrid::GammaTubeIterator MyOpenGammaTubeIterator;
typedef MyOpenGrid::ZeroCrossingIterator MyOpenZeroCrossingIterator;
typedef MyOpenGrid::InitParams MyOpenInitParams;

typedef LevelSet::LevelSet3D<MyOpenGrid> MyOpenLevelSet;
typedef MyOpenLevelSet::InitParams MyOpenLSInitParams;
typedef MyOpenLevelSet::SimulationParameters MyOpenSimulationParameters;
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
typedef DTGridTraitsClosedLevelSet<Data, Index, Real, UInt> MyTraitsClosedLevelSet;

typedef Grids::DTGrid<MyTraitsClosedLevelSet> MyGrid;
typedef MyGrid::TubeIterator MyTubeIterator;
typedef MyGrid::EntireTubeIterator MyEntireTubeIterator;
typedef MyGrid::BetaTubeIterator MyBetaTubeIterator;
typedef MyGrid::GammaTubeIterator MyGammaTubeIterator;
typedef MyGrid::ZeroCrossingIterator MyZeroCrossingIterator;
typedef MyGrid::InitParams MyInitParams;

typedef LevelSet::LevelSet3D<MyGrid> MyLevelSet;
typedef MyLevelSet::InitParams MyLSInitParams;
typedef MyLevelSet::SimulationParameters MySimulationParameters;
////////////////////////////////////////////////////////////////////////////////////////


void meanValue();
void augmentedLevelSet();
void randomAccess();
void averageMeanCurvature();
void localizedLevelSetComputation();
void randomAccessOpenLevelSet();
void rebuildOpenLevelSet();
void csgOperation();
void csgOperationWithTransform();
void csgOperationOpenLevelSet();
void neighborAccessOpenLevelSet();
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
    int testType;

    if (argc != 3)
    {
	cout << "Usage: DTGrid_Test <svol-file-name> <test-id>" << endl;
	cout << "Tests: " << endl;
	cout << " 1: Iteration (compute mean value)" << endl;
	cout << " 2: Stencil Iteration (compute average mean curvature)" << endl;
	cout << " 3: Random Access (compute mean value)" << endl;
	cout << " 4: Augmented Level Set Computations" << endl; 
	cout << " 5: Localized Level Set Computations on Open Level Sets" << endl;
	cout << " 6: Random Access in Open Level Sets" << endl;
	cout << " 7: Narrow Band Rebuild in Open Level Sets" << endl;
	cout << " 8: CSG Operations" << endl;
	cout << " 9: CSG Operations that take Transformations into Account" << endl;
	cout << "10: CSG Operations on Open Level Sets" << endl;
	cout << "11: Neighbor Access in Open Level Sets" << endl;
	cout << "12: Closest Point" << endl;
	cout << "13: Perform Unit Tests" << endl;
	exit(-1);
    }

    inputFileName = argv[1];
    testType = atoi(argv[2]);

    switch (testType)
    {
    case 1:
	meanValue();
	break;
    case 2:
	averageMeanCurvature();
	break;
    case 3:
	randomAccess();
	break;
    case 4:
	augmentedLevelSet();
	break;
    case 5:
	localizedLevelSetComputation();
	break;
    case 6:
	randomAccessOpenLevelSet();
	break;
    case 7:
	rebuildOpenLevelSet();
	break;
    case 8:
	csgOperation();
	break;
    case 9:
	csgOperationWithTransform();
	break;
    case 10:
	csgOperationOpenLevelSet();
	break;
    case 11:
	neighborAccessOpenLevelSet();
	break;
    case 12:
	closestPoint();
	break;
    case 13:
	meanValue();
	averageMeanCurvature();
	randomAccess();
	augmentedLevelSet();
	localizedLevelSetComputation();
	randomAccessOpenLevelSet();
	rebuildOpenLevelSet();
	csgOperation();
	csgOperationWithTransform();
	csgOperationOpenLevelSet();
	neighborAccessOpenLevelSet();
	closestPoint();
	break;
    }

    return 1;
}




void rebuildOpenLevelSet()
{
    MyOpenGrid *grid, *subset;
    MyOpenInitParams initParams;
    Index bbox[3][2];
    Index bboxDim[3];
    UInt numValuesInOriginalSubset;
    UInt numValuesInRebuiltSubset;


    grid = new MyOpenGrid(initParams);
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
    subset = grid->copy(bbox);

    subset = grid->copy(bbox);
    // set the open level set bbox on the subset.
    // this ensures that the level set does not expand beyond
    // this bounding box during dilation and rebuild operations
    subset->setOpenLevelSetBBox(bbox);

    // save subset as sparse volume for inspection
    subset->saveSparseVolume("rebuildopenlevelset_subset.svol");

    numValuesInOriginalSubset = subset->getNumValues();
    cout << "Number of values in original subset = " << numValuesInOriginalSubset << endl;
    // The assumption of the rebuildNarrowBandMethod is that the narrow band was fully contained 
    // within the openLevelSetBoundingBox prior to the rebuild 
    subset->rebuildNarrowBand();
    numValuesInRebuiltSubset = subset->getNumValues();
    cout << "Number of values in rebuilt subset = " << numValuesInRebuiltSubset << endl;

    subset->saveSparseVolume("rebuildopenlevelset_subset_rebuilt.svol");

    if (numValuesInRebuiltSubset < numValuesInOriginalSubset)
    {
	cout << "ERROR in rebuildOpenLevelSet():" << endl;
	cout << " Number of values in original subset <  Number of values in rebuilt subset" << endl;
    }
    else
    {
	cout << "No runtime errors in rebuildOpenLevelSet(), check result!" << endl;
    }

    delete grid;
    delete subset;
}




void randomAccessOpenLevelSet()
{
    MyOpenGrid *grid, *subset;
    MyOpenInitParams initParams;
    Index bbox[3][2];
    Index bboxDim[3];
    Index x, y, z;
    bool error = false;

    grid = new MyOpenGrid(initParams);
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
		Real value1 = (*grid)(x,y,z);
		Real value2 = (*subset)(x,y,z);  // random access
		Real value3;
		(*subset)(x,y,z,&value3);  // random access, approach 2

		if (value1 != value2 || value1 != value3)
		{
		    cout << "ERROR in randomAccessOpenLevelSet() at (" << x << "," << y << "," << z << ")" << endl;
		    cout << "value1 = " << value1 << endl;
		    cout << "value2 = " << value2 << endl;
		    cout << "value3 = " << value3 << endl;
		    error = true;
		}
	    }
	}
    }
    if (!error)
    {
	cout << "No errors in randomAccessOpenLevelSet()!" << endl;
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
void localizedLevelSetComputation()
{
    typedef LevelSet::MeanCurvatureFlowScalarField3D<Real, Index, MyOpenGrid, MyOpenGammaTubeIterator> MyScalarField;

    MyOpenGrid *grid, *subset;
    MyOpenInitParams initParams;
    Index bbox[3][2];
    MyOpenLSInitParams lsInitParams;
    MyOpenLevelSet *ls;
    MyOpenSimulationParameters simParameters;
    MyScalarField *sf;


    grid = new MyOpenGrid(initParams);
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
    ls->propagateDt<Grids::SF_FIRSTORDER_CURVATURE, Grids::SF_FIRSTORDER, MyScalarField >(sf, simParameters);	

    delete ls;
    delete sf;

    // paste the subset narrow band back into the original level set
    grid->paste(subset, bbox);
    delete subset;

    // save the final result as a sparse volume (svol) file
    grid->saveSparseVolume("localizedLevelSetComputation_final_result.svol");
    delete grid;

    cout << "No runtime errors in localizedLevelSetComputation(), check result!" << endl;
}






/**
 * This example demonstrates iterating over the narrow band to compute the mean value 
 */
void meanValue()
{
    MyGrid *grid;
    MyInitParams initParams;
    MyTubeIterator iter, iend;
    Real mean = 0;

    grid = new MyGrid(initParams);
    grid->loadSparseVolume(inputFileName);

    iter = grid->beginTubeIterator();
    iend = grid->endTubeIterator();

    while (iter != iend)
    {
	Real value = iter.getValue();
	mean += value;
	iter++;
    }

    mean /= grid->getNumValues();

    cout << "mean = " << mean << endl;

    cout << "No runtime errors in meanValue(), check result!" << endl;
}


/** 
 * This example demonstrates the random access operator to compute the mean value
 * and compares the result to the mean value computed using iteration
 */ 
void randomAccess()
{
    MyGrid *grid;
    MyInitParams initParams;
    MyTubeIterator iter, iend;
    Real value, value2;
    Real meanIter = 0;
    Real meanRandom = 0;
    Index x=0, y=0, z=0;
    bool error = false;

    grid = new MyGrid(initParams);
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
	cout << "ERROR in randomAccess(), value = " << value << ", value2 = " << value2 << " at (" << x << "," << y << "," << z << ")" << endl;
    }
    else
    {
	if (meanIter != meanRandom)
	{
	    cout << "ERROR in randomAccess(), meanIter != meanRandom" << endl;
	    cout << "meanIter = " << meanIter << endl;
	    cout << "meanRandom = " << meanRandom << endl;
	}
	else
	{
	    cout << "No error in randomAccess()!" << endl;
	}
    }

    meanIter /= grid->getNumValues();

    cout << "meanIter = " << meanIter << endl;
}



/**
 *  Augmenting a level set with additional data is useful in many applications.
 *  Augmented data can be e.g. texture, colors, normals or velocities from a fluid simulation.
 *  This example shows how the DT-Grid can be augmented with additional data arrays and how to maintain
 *  this data during narrow band rebuild.
 *  In this simple example, the data is simple an array of 0 on the surface, +1 outside the interface and -1 inside the interface. */
void augmentedLevelSet()
{
    MyGrid *grid;
    MyInitParams initParams;
    UInt size, sizeAfterRebuild;
    UInt *permutationArray;
    int *augmentedData, *augmentedDataAfterRebuild;
    bool error = false;

    grid = new MyGrid(initParams);
    grid->loadSparseVolume(inputFileName);

    size = grid->getNumValues();
    augmentedData = new int[size+2];  // make room for inside/outside values
    for (UInt i=0; i<size; i++)
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
    for (UInt i=0; i<sizeAfterRebuild; i++)
    {
	augmentedDataAfterRebuild[i] = augmentedData[ permutationArray[i] ];
    }
    // set inside value
    augmentedDataAfterRebuild[sizeAfterRebuild] = -1;
    // set outside value
    augmentedDataAfterRebuild[sizeAfterRebuild+1] = 1;	

    // Test
    for (UInt i=0; i<sizeAfterRebuild && !error; i++)
    {
	if (permutationArray[i]>size+2)
	{
	    cout << "ERROR in augmentedLevelSet() test 1 at i = " << i << endl;
	    cout << "permutation array value = " << permutationArray[i] << endl;
	    error = true;
	}
    }
    // Test
    for (UInt i=0; i<sizeAfterRebuild && !error; i++)
    {
	if (grid->getValue(i) < 0 && augmentedDataAfterRebuild[i]>=0 ||
	    grid->getValue(i) > 0 && augmentedDataAfterRebuild[i]<=0 ||
	    grid->getValue(i) == 0 && augmentedDataAfterRebuild[i]!=0)
	{
	    cout << "ERROR in augmentedLevelSet() test 2 at i = " << i << endl;
	    cout << "Grid value = " << grid->getValue(i) << endl;
	    cout << "Augmented data value = " << augmentedDataAfterRebuild[i] << endl;
	    cout << "sizeAfterRebuild = " << sizeAfterRebuild << endl;
	    cout << "permutation array value = " << permutationArray[i] << endl;
	    error = true;
	}
    }
    if (!error)
    {
	cout << "No errors in augmentedLevelSet()!" << endl;
    }


    delete[] augmentedData;
    delete[] permutationArray;
    delete grid;
    delete[] augmentedDataAfterRebuild;
}






/** This example demonstrates iteration with a stencil over the narrow band 
 *  to compute the average mean curvature */
void averageMeanCurvature()
{
    MyGrid *grid;
    MyEntireTubeIterator iter, iend;
    MyInitParams initParams;

    grid = new MyGrid(initParams);
    grid->loadSparseVolume(inputFileName);

    // compute average mean curvature
    Real dx = grid->getScale();
    Real epsilon = (Real)(1.5 * dx);
    Real volElement = dx * dx * dx;
    Real deltaVal;
    Real c1 = (Real)(1.0 / (Real(2.0) * epsilon));
    Real c2 = (Real)(M_PI / epsilon);
    Real H;
    Real gradLen;
    Real area;
    Real cellArea;
    Real phiVal;
    Real averageMeanCurvature;

    iter = grid->beginEntireTube<Grids::SF_FIRSTORDER_CURVATURE, false>();
    iend = grid->endEntireTube<Grids::SF_FIRSTORDER_CURVATURE, false>();

    averageMeanCurvature = 0;
    area = 0;

    while (iter != iend)
    {
	phiVal = iter.getValue();

	if (fabs(phiVal) <= epsilon)
	{
	    deltaVal = c1 + c1 * cos( phiVal * c2 );
	    H = iter.meanCurvature<Grids::SF_FIRSTORDER_CURVATURE>();
	    gradLen = iter.gradientLength();
	    cellArea = deltaVal * gradLen * volElement;
	    area += cellArea;
	    averageMeanCurvature +=  H * cellArea; 
	}

	iter.operator++<Grids::SF_FIRSTORDER_CURVATURE, false>();
    }

    averageMeanCurvature /= area;

    cout << "Average mean Curvature is " << averageMeanCurvature << endl;

    cout << "No runtime errors in averageMeanCurvature(), check result!" << endl;

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
void csgOperation()
{
    MyGrid *grid1, *grid2;
    MyEntireTubeIterator iter, iend;
    MyInitParams initParams;
    Index bboxDim[3];

    grid1 = new MyGrid(initParams);
    grid1->loadSparseVolume(inputFileName);
    grid1->boundingBoxDim(bboxDim);
    bboxDim[0] /= 8;
    bboxDim[1] /= 8;
    bboxDim[2] /= 8;
    grid1->translateGrid(bboxDim);

    // Union
    grid2 = new MyGrid(initParams);
    grid2->loadSparseVolume(inputFileName);
    cout << "Union.."; cout.flush();
    grid2->CSG_GA(*grid1, MyGrid::CSG_UNION);
    cout << "done" << endl;
    grid2->saveSparseVolume("union_result.svol");
    delete grid2;

    // Intersection
    grid2 = new MyGrid(initParams);
    grid2->loadSparseVolume(inputFileName);
    cout << "Intersection.."; cout.flush();
    grid2->CSG_GA(*grid1, MyGrid::CSG_INTERSECTION);
    cout << "done" << endl;
    grid2->saveSparseVolume("intersection_result.svol");
    delete grid2;

    // Difference
    grid2 = new MyGrid(initParams);
    grid2->loadSparseVolume(inputFileName);
    cout << "Difference.."; cout.flush();
    grid2->CSG_GA(*grid1, MyGrid::CSG_DIFFERENCE);
    cout << "done" << endl;
    grid2->saveSparseVolume("difference_result.svol");
    delete grid2;

    delete grid1;

    cout << "No runtime errors in csgOperation(), check result!" << endl;
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
void csgOperationWithTransform()
{
    MyGrid *grid1, *grid2;
    MyEntireTubeIterator iter, iend;
    MyInitParams initParams;

    Math::TrilinearInterpolator<Real, Real> triInterp;

    grid1 = new MyGrid(initParams);
    grid1->loadSparseVolume(inputFileName);
    grid1->setRotationXYZ(1,1,1);  // angles specified in radians

    // Union
    grid2 = new MyGrid(initParams);
    grid2->loadSparseVolume(inputFileName);
    cout << "Union.."; cout.flush();
    grid2->CSG(*grid1, MyGrid::CSG_UNION, &triInterp);
    cout << "done" << endl;
    grid2->saveSparseVolume("union_transform_result.svol");
    delete grid2;

    // Intersection
    grid2 = new MyGrid(initParams);
    grid2->loadSparseVolume(inputFileName);
    cout << "Intersection.."; cout.flush();
    grid2->CSG(*grid1, MyGrid::CSG_INTERSECTION, &triInterp);
    cout << "done" << endl;
    grid2->saveSparseVolume("intersection_transform_result.svol");
    delete grid2;

    // Difference
    grid2 = new MyGrid(initParams);
    grid2->loadSparseVolume(inputFileName);
    cout << "Difference.."; cout.flush();
    grid2->CSG(*grid1, MyGrid::CSG_DIFFERENCE, &triInterp);
    cout << "done" << endl;
    grid2->saveSparseVolume("difference_transform_result.svol");
    delete grid2;

    delete grid1;

    cout << "No runtime errors in csgOperationWithTransform(), check result!" << endl;
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
void csgOperationOpenLevelSet()
{
    MyOpenGrid *grid1, *grid2, *subset1, *subset2, *subset3;
    MyOpenTubeIterator iter, iend, iter2, iend2;
    MyOpenInitParams initParams;
    Index bboxDim[3];
    Index bbox[3][2];

    grid1 = new MyOpenGrid(initParams);
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
    grid2 = new MyOpenGrid(initParams);
    grid2->loadSparseVolume(inputFileName);
    subset2 = grid2->copy(bbox);
    subset2->setOpenLevelSetBBox(bbox);
    cout << "Union.."; cout.flush();
    subset2->saveSparseVolume("subset2.svol");
    subset2->CSG_GA(*subset1, MyOpenGrid::CSG_UNION);
    cout << "done" << endl;
    subset2->saveSparseVolume("union_openls_result.svol");

    // Intersection
    grid2 = new MyOpenGrid(initParams);
    grid2->loadSparseVolume(inputFileName);
    subset2 = grid2->copy(bbox);
    subset2->setOpenLevelSetBBox(bbox);
    cout << "Intersection.."; cout.flush();
    subset2->CSG_GA(*subset1, MyOpenGrid::CSG_INTERSECTION);
    cout << "done" << endl;
    subset2->saveSparseVolume("intersection_openls_result.svol");
    delete grid2;
    delete subset2;

    // Difference
    grid2 = new MyOpenGrid(initParams);
    grid2->loadSparseVolume(inputFileName);
    subset2 = grid2->copy(bbox);
    subset2->setOpenLevelSetBBox(bbox);
    cout << "Difference.."; cout.flush();
    subset2->CSG_GA(*subset1, MyOpenGrid::CSG_DIFFERENCE);
    cout << "done" << endl;
    subset2->saveSparseVolume("difference_openls_result.svol");


    // Compare to difference with closed level set followed by a copy operation.
    // In particular this should give the same result.
    cout << "Difference with closed level set.."; cout.flush();
    grid2->CSG_GA(*grid1, MyOpenGrid::CSG_DIFFERENCE);
    cout << "done" << endl;
    subset3 = grid2->copy(bbox);

    iter = subset2->beginTubeIterator();
    iend = subset2->endTubeIterator();
    iter2 = subset3->beginTubeIterator();
    iend2 = subset3->endTubeIterator();

    int numPointsOk = 0;
    bool error = false;

    while ((iter!=iend || iter2!=iend2) && !error)
    {
	Index c[3], c2[3];

	iter.getIndex(&c[0], &c[1], &c[2]);
	iter2.getIndex(&c2[0], &c2[1], &c2[2]);

	if (c[0] != c2[0] || c[1] != c2[1] || c[2] != c2[2])
	{
	    error = true;
	    cout << "ERROR in csgOperationOpenLevelSet():" << endl;
	    cout << "point from open level set: " << c[0] << " " << c[1] << " " << c[2] << endl;
	    cout << "point from closed level set: " << c2[0] << " " << c2[1] << " " << c2[2] << endl;
	    cout << "NumPointsOK = " << numPointsOk << endl;
	    cout << endl;
	}
	else
	{
	    numPointsOk++;
	}

	if (iter!=iend)
	    iter++;
	if (iter2!=iend2)
	    iter2++;
    }
    if (!error)
    {
	cout << "No errors in csgOperationOpenLevelSet()!" << endl;
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
void neighborAccessOpenLevelSet()
{
    MyOpenGrid *grid, *subset;
    MyOpenInitParams initParams;
    Index bbox[3][2];
    Index bboxDim[3];
    Index x, y, z;
    bool error = false;

    grid = new MyOpenGrid(initParams);
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
		    MyOpenGrid::Locator loc;

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
			cout << "ERROR in neighborAccessOpenLevelSet() at (" << x << "," << y << "," << z << ")" << endl;
			cout << "xm1 = " << xm << endl;
			cout << "xm2 = " << xm2 << endl;
			error = true;
		    }
		    if (xp!=xp2)
		    {
			cout << "ERROR in neighborAccessOpenLevelSet() at (" << x << "," << y << "," << z << ")" << endl;
			cout << "xp1 = " << xp << endl;
			cout << "xp2 = " << xp2 << endl;
			error = true;
		    }

		    if (ym!=ym2)
		    {
			cout << "ERROR in neighborAccessOpenLevelSet() at (" << x << "," << y << "," << z << ")" << endl;
			cout << "ym1 = " << ym << endl;
			cout << "ym2 = " << ym2 << endl;
			error = true;
		    }
		    if (yp!=yp2)
		    {
			cout << "ERROR in neighborAccessOpenLevelSet() at (" << x << "," << y << "," << z << ")" << endl;
			cout << "yp1 = " << yp << endl;
			cout << "yp2 = " << yp2 << endl;
			error = true;
		    }

		    if (zm!=zm2)
		    {
			cout << "ERROR in neighborAccessOpenLevelSet() at (" << x << "," << y << "," << z << ")" << endl;
			cout << "zm1 = " << zm << endl;
			cout << "zm2 = " << zm2 << endl;
			error = true;
		    }
		    if (zp!=zp2)
		    {
			cout << "ERROR in neighborAccessOpenLevelSet() at (" << x << "," << y << "," << z << ")" << endl;
			cout << "zp1 = " << zp << endl;
			cout << "zp2 = " << zp2 << endl;
			error = true;
		    }

		}
	    }
	}
    }
    if (!error)
    {
	cout << "No errors in neighborAccessOpenLevelSet()!" << endl;
    }

    delete grid;
    delete subset;
}




void closestPoint()
{
    MyGrid *grid;
    MyInitParams initParams;
    Index bbox[3][2];
    Index x, y, z;
    Index p[3];               // grid point
    Real n[3];               // normal vector at grid point
    Real cd;                  // closest distance to interface at grid point
    bool inside;              // true if grid point is inside the object
    bool inNarrowBand;        // true if grid point is inside the narrow band

    UInt numInNarrowBand = 0;
    UInt numInsideObject = 0;

    grid = new MyGrid(initParams);
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

    cout << "numInsideObject = " << numInsideObject << endl;
    cout << "numInNarrowBand = " << numInNarrowBand << endl;

    if (numInNarrowBand != grid->getNumValues())
    {
	cout << "ERROR in closestPoint(), numInNarrowBand = " << numInNarrowBand << ", grid->getNumValues() = " << grid->getNumValues() << endl;
    }
    else
    {
	cout << "No errors in closestPoint()!" << endl;
    }

}
