#!/bin/bash

# Build script for the integrated Dense-PCE + EBBkC system
# Builds EBBkC as a library and links it with a specified Dense-PCE source file.
#
# Usage:
#   ./build_integrated.sh <source_file.cpp>          # release (default)
#   ./build_integrated.sh <source_file.cpp> release  # explicit release
#   ./build_integrated.sh <source_file.cpp> debug    # debug + CUI invariant on
#
# Release build:
#   -O3 -march=native -DNDEBUG
#   DENSE_PCE_CHECK_DPSI expands to 0, invariant is compiled out.
#   Use this for all performance measurements.
#
# Debug build:
#   -O0 -g -DDENSE_PCE_CHECK_DPSI=1
#   Invariant runs on every iter() call in CUI mode.
#   Use this only to verify correctness on small graphs.

set -e  # Exit on any error

# Check for command line argument
if [ $# -eq 0 ]; then
    echo "Usage: $0 <source_file.cpp> [release|debug]"
    echo "Example: $0 dense-pce-cui.cpp           # release (default)"
    echo "Example: $0 dense-pce-cui.cpp debug     # debug with CUI invariant"
    exit 1
fi

SOURCE_FILE="$1"
BUILD_MODE="${2:-release}"

# Check if source file exists
if [ ! -f "$SOURCE_FILE" ]; then
    echo "Error: Source file '$SOURCE_FILE' not found"
    exit 1
fi

# Select compile flags per mode
case "$BUILD_MODE" in
    release)
        DENSE_PCE_CXXFLAGS="-O3 -march=native -DNDEBUG"
        SUFFIX=""
        ;;
    debug)
        DENSE_PCE_CXXFLAGS="-O0 -g -DDENSE_PCE_CHECK_DPSI=1"
        SUFFIX="-debug"
        ;;
    *)
        echo "Error: unknown build mode '$BUILD_MODE' (expected 'release' or 'debug')"
        exit 1
        ;;
esac

# Extract base name for output executable
BASE_NAME=$(basename "$SOURCE_FILE" .cpp)
OUTPUT_NAME="${BASE_NAME}-integrated${SUFFIX}"

echo "=== Building Integrated Dense-PCE + EBBkC System ==="
echo "Source file : $SOURCE_FILE"
echo "Build mode  : $BUILD_MODE"
echo "Flags       : $DENSE_PCE_CXXFLAGS"
echo "Executable  : $OUTPUT_NAME"

# Clean up old build artifacts to avoid path conflicts
echo "=== Cleaning up old build artifacts ==="
rm -rf EBBkC/src/build

# Create build directory if it doesn't exist
mkdir -p build_integrated
cd build_integrated

echo "=== Step 1: Building EBBkC as a library ==="

# Build EBBkC core library (always release — EBBkC isn't what we're debugging)
cd ../EBBkC/src
mkdir -p build
cd build

# Clean any existing CMake cache to avoid path conflicts
rm -f CMakeCache.txt
rm -rf CMakeFiles

echo "Configuring EBBkC with CMake..."
cmake .. -DCMAKE_BUILD_TYPE=Release
echo "Building EBBkC..."
make -j$(nproc)

# Copy the library to main directory
cp libebbkc_core.a ../../../build_integrated/
cd ../../../build_integrated

echo "=== Step 2: Compiling $SOURCE_FILE with EBBkC library ==="

# Compile with EBBkC. Two attempts for different OpenMP library names.
echo "Compiling integrated executable..."
g++ -std=c++17 $DENSE_PCE_CXXFLAGS -fopenmp \
    -I../EBBkC/src \
    -I../EBBkC/src/truss/dependencies/sparsepp \
    -I../EBBkC/src/truss/dependencies/libpopcnt \
    -I../EBBkC/src/truss/dependencies \
    -I../EBBkC/src/truss/util \
    -I../EBBkC/src/truss/decompose \
    ../$SOURCE_FILE \
    libebbkc_core.a \
    ../EBBkC/src/build/libcommon-utils.a \
    ../EBBkC/src/build/libgraph-pre-processing.a \
    -ltbb -lpthread \
    -o $OUTPUT_NAME 2>/dev/null || \
g++ -std=c++17 $DENSE_PCE_CXXFLAGS -fopenmp \
    -I../EBBkC/src \
    -I../EBBkC/src/truss/dependencies/sparsepp \
    -I../EBBkC/src/truss/dependencies/libpopcnt \
    -I../EBBkC/src/truss/dependencies \
    -I../EBBkC/src/truss/util \
    -I../EBBkC/src/truss/decompose \
    ../$SOURCE_FILE \
    libebbkc_core.a \
    ../EBBkC/src/build/libcommon-utils.a \
    ../EBBkC/src/build/libgraph-pre-processing.a \
    -lgomp -ltbb -lpthread \
    -o $OUTPUT_NAME

echo "=== Build completed successfully! ==="
echo "Executable: ./build_integrated/$OUTPUT_NAME"
echo ""
echo "Quick sanity check — invariant should be OFF in release builds:"
echo "  strings build_integrated/$OUTPUT_NAME | grep -c 'D_PSI invariant failed'"
echo "  (expected: 0 for release, 1 for debug)"
echo ""
echo "Usage example:"
echo "  ./build_integrated/$OUTPUT_NAME testGraph/gplus/gplus.grh --minimum 10 --theta 0.9 --mode 5"