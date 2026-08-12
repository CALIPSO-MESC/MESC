#!/bin/bash
# Shell script for building MESC.

set -e

# --------------------------------------------------
# Load environment
#
# This will be specific to the machine you are working on and will need to be
# completed by the user. For example:
# ```
# module purge
# module load oneapi23u1   # load intel compiler
# module load netcdf_intel # load netcdf library
# ```
# --------------------------------------------------

# --------------------------------------------------
# Setup build directory
# --------------------------------------------------

BUILD_DIR=build
if [ ! -d "${BUILD_DIR}" ]; then
  mkdir "${BUILD_DIR}"
fi
cd "${BUILD_DIR}"

# --------------------------------------------------
# Configure & build
# --------------------------------------------------

cmake .. -DCMAKE_Fortran_COMPILER=ifort
cmake --build . -j

echo "Build successful."
echo "The $(main) executable is available in both test/ and build/"
