#!/bin/bash

# Build script for the integrated Dense-PCE + EBBkC system
# This script compiles EBBkC as a library and links it with dense-pce-mod-edge-order

set -e  # Exit on any error

echo "=== Building Integrated Dense-PCE + EBBkC System ==="

# Check if we're in the right directory
if [ ! -f "dense-pce-mod-edge-order.cpp" ]; then
    echo "Error: Please run this script from the Dense-PCE-main directory"
    exit 1
fi

# Clean up old build artifacts to avoid path conflicts
echo "=== Cleaning up old build artifacts ==="
rm -rf build_integrated
rm -rf EBBkC/src/build

# Create build directory
mkdir -p build_integrated
cd build_integrated

echo "=== Step 1: Building EBBkC as a library ==="

# Build EBBkC core library
cd ../EBBkC/src
mkdir -p build
cd build

# Clean any existing CMake cache to avoid path conflicts
rm -f CMakeCache.txt
rm -rf CMakeFiles

# Configure and build EBBkC
echo "Configuring EBBkC with CMake..."
cmake .. -DCMAKE_BUILD_TYPE=Release
echo "Building EBBkC..."
make -j$(nproc)

# Copy the library to main directory
cp libebbkc_core.a ../../../build_integrated/
cd ../../../build_integrated

echo "=== Step 2: Compiling dense-pce-mod-edge-order with EBBkC library ==="

# Compile dense-pce-mod-edge-order with EBBkC library
# Try different OpenMP library names for different systems
echo "Compiling integrated executable..."
g++ -std=c++17 -O3 -march=native -fopenmp \
    -I../EBBkC/src \
    -I../EBBkC/src/truss/dependencies/sparsepp \
    -I../EBBkC/src/truss/dependencies/libpopcnt \
    -I../EBBkC/src/truss/dependencies \
    -I../EBBkC/src/truss/util \
    -I../EBBkC/src/truss/decompose \
    ../dense-pce-mod-edge-order.cpp \
    libebbkc_core.a \
    ../EBBkC/src/build/libcommon-utils.a \
    ../EBBkC/src/build/libgraph-pre-processing.a \
    -ltbb -lpthread \
    -o dense-pce-mod-edge-order-integrated 2>/dev/null || \
g++ -std=c++17 -O3 -march=native -fopenmp \
    -I../EBBkC/src \
    -I../EBBkC/src/truss/dependencies/sparsepp \
    -I../EBBkC/src/truss/dependencies/libpopcnt \
    -I../EBBkC/src/truss/dependencies \
    -I../EBBkC/src/truss/util \
    -I../EBBkC/src/truss/decompose \
    ../dense-pce-mod-edge-order.cpp \
    libebbkc_core.a \
    ../EBBkC/src/build/libcommon-utils.a \
    ../EBBkC/src/build/libgraph-pre-processing.a \
    -lgomp -ltbb -lpthread \
    -o dense-pce-mod-edge-order-integrated

echo "=== Build completed successfully! ==="
echo "Executable: ./build_integrated/dense-pce-mod-edge-order-integrated"
echo ""
echo "Usage example:"
echo "./build_integrated/dense-pce-mod-edge-order-integrated testGraphs/gplus/gplus.grh --minimum 10 --theta 0.9"
