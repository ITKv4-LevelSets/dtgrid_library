project( dtgrid )

cmake_minimum_required(VERSION 2.8)

include_directories( BEFORE
  ${dtgrid_SOURCE_DIR}/Core/Exception
  ${dtgrid_SOURCE_DIR}/Algorithms/LevelSet
  ${dtgrid_SOURCE_DIR}/Algorithms/LevelSet/ScalarFields
  ${dtgrid_SOURCE_DIR}/Algorithms/LevelSet/VelocityFields
  ${dtgrid_SOURCE_DIR}/Algorithms/Math
  ${dtgrid_SOURCE_DIR}/DataStructures/Arrays
  ${dtgrid_SOURCE_DIR}/DataStructures/Containers
  ${dtgrid_SOURCE_DIR}/DataStructures/Grids
  ${dtgrid_SOURCE_DIR}/DataStructures/Grids/IO
  ${dtgrid_SOURCE_DIR}/DataStructures/Grids/DTGrid
  ${dtgrid_SOURCE_DIR}/DataStructures/Matrices
  ${dtgrid_SOURCE_DIR}/DataStructures/Vectors
)

add_executable( test Test/DTGrid_Test/DTGrid_Test.cpp )

