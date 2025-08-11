# Dense-PCE Project Structure

## Overview
This repository contains implementations of various pseudo-clique enumeration algorithms, including Dense-PCE, PCE, FPCE, and EBBkC. The project focuses on finding dense subgraphs (pseudo-cliques) in large-scale networks with different optimization strategies.

## Repository Structure

```
CSE491-/
├── README.md                                    # Main project README
├── WSL-Ubuntu-cmds.txt                         # WSL Ubuntu commands reference
└── Dense_PCE/                                  # Main project directory
    └── Dense-PCE-main/                         # Core implementation directory
        │
        ├── Core Implementations/               # Main algorithm implementations
        │   ├── dense-pce.cpp                   # Original Dense-PCE implementation (C++)
        │   ├── dense-pce-mod.cpp               # Modified Dense-PCE with Turan filtering & maximal enumeration
        │   ├── dense-pce                       # Compiled Dense-PCE executable
        │   ├── dense-pce-mod                   # Compiled modified Dense-PCE executable
        │   ├── pce                             # Compiled PCE executable
        │   └── BBkC                            # Compiled EBBkC executable
        │
        ├── Algorithm Source Directories/       # Algorithm source code
        │   ├── pce12/                          # Original PCE implementation (C)
        │   │   ├── pce.c                       # Main PCE algorithm with maximal enumeration
        │   │   ├── problem.c/h                 # Problem definition and parameters
        │   │   ├── sgraph.c/h                  # Graph data structures
        │   │   ├── queue.c/h                   # Queue implementation
        │   │   ├── itemset.c/h                 # Itemset handling and output
        │   │   ├── vec.c/h                     # Vector operations
        │   │   ├── aheap.c/h                   # Heap implementation
        │   │   ├── alist.c/h                   # Adjacency list management
        │   │   ├── base.c/h                    # Base utilities
        │   │   ├── undo.c/h                    # Undo functionality
        │   │   ├── stdlib2.c/h                 # Extended standard library
        │   │   ├── makefile                    # Build configuration
        │   │   ├── readme.txt                  # PCE documentation
        │   │   ├── test*.grh                   # Test graph files
        │   │   └── *.pl                        # Perl utilities for graph processing
        │   │
        │   ├── FPCE/                           # Fast Pseudo-Clique Enumeration
        │   │   ├── code/                       # FPCE source code
        │   │   │   ├── FPCE/                   # Main FPCE implementation (33 files)
        │   │   │   ├── PCE/                    # PCE variant (23 files)
        │   │   │   ├── ODES/                   # Dense subgraph algorithm (8 files)
        │   │   │   ├── SIGMOD24-MQCE-main/     # FastQC implementation (24 files)
        │   │   │   └── cluster_one-1.0.jar     # ClusterONE tool
        │   │   ├── real-synth-graph-analysis/  # Graph analysis experiments
        │   │   │   ├── collect_graph_info.ipynb
        │   │   │   ├── command_generator.ipynb
        │   │   │   ├── generate_synth_graphs.ipynb
        │   │   │   ├── real-graphs/            # Real-world graph datasets
        │   │   │   └── synth-graphs/           # Synthetic graph information
        │   │   ├── e-coli-analysis/            # E. coli biological network analysis
        │   │   │   ├── e-coli-data/            # E. coli protein interaction data
        │   │   │   ├── evaluate_ecoli_pseudocliques.ipynb
        │   │   │   ├── evaluate_ecoli_pseudocliques_vs_fastqc.ipynb
        │   │   │   ├── process_ecoli_data.sh
        │   │   │   └── subgraph_densities.py
        │   │   ├── README.md                   # FPCE documentation
        │   │   ├── Graph_check.py              # Graph preprocessing utility
        │   │   └── *.pl                        # Graph conversion utilities
        │   │
        │   └── EBBkC/                          # Efficient k-clique listing with GPU support
        │       ├── src/                        # EBBkC source code
        │       │   ├── main.cpp                # Main EBBkC entry point
        │       │   ├── edge_oriented.cpp/h     # Edge-oriented branching algorithm
        │       │   ├── set_operation.cpp/h     # Set operations
        │       │   ├── def.cpp/h               # Definitions and configurations
        │       │   ├── CMakeLists.txt          # CMake build configuration
        │       │   ├── build/                  # Build directory (100 files)
        │       │   └── truss/                  # Truss decomposition (344 files)
        │       │       ├── decompose/          # Truss decomposition algorithms
        │       │       ├── dependencies/       # External libraries
        │       │       │   ├── cub/            # CUDA utilities
        │       │       │   ├── libpopcnt/      # Population count library
        │       │       │   ├── moderngpu/      # Modern GPU utilities
        │       │       │   ├── sparsehash-yche/# Sparse hash implementation
        │       │       │   └── sparsepp/       # Sparse++ library
        │       │       └── util/               # Graph utilities and reordering
        │       ├── dataset/                    # Test datasets
        │       │   ├── facebook/               # Facebook network data
        │       │   ├── facebook_100/           # Facebook subset (100 nodes)
        │       │   ├── nasasrb/                # NASA collaboration network
        │       │   ├── test10/test15/test20/   # Small test graphs
        │       │   ├── testK4/testK7/          # K-clique test datasets
        │       │   └── *.clean, *.edges        # Processed graph files
        │       ├── Cohesive_subgraph_book/     # Graph algorithms reference implementation
        │       │   ├── core_decomposition/     # Core decomposition algorithms
        │       │   ├── data_structures/        # Graph data structures
        │       │   ├── datasets/               # Example datasets
        │       │   ├── densest_subgraph/       # Densest subgraph algorithms
        │       │   ├── edge_connectivity_decomposition/ # Edge connectivity algorithms
        │       │   └── utilities/              # Utility headers
        │       ├── README.md                   # EBBkC documentation
        │       ├── EBBkC_TR.pdf                # Technical report
        │       ├── make_100.py                 # Graph generation utility
        │       ├── test20.py                   # Testing script
        │       └── WSL-Ubuntu-commands.txt     # Ubuntu build instructions
        │
        ├── Graph Datasets/                     # Graph data files
        │   ├── testGraphs/                     # Test graph files
        │   │   └── fpce_graph/                 # FPCE test graphs with K-clique files
        │   │       ├── fpce_graph.grh          # Main test graph
        │   │       ├── fpce_graph.txt/.clean/.edges # Various formats
        │   │       ├── cliques_K1/K2/K3/K4     # Pre-computed clique files
        │   │       └── b_adj.bin, b_degree.bin # Binary graph representation
        │   ├── synth-graphs-1000/              # Synthetic scale-free graphs (1000 nodes)
        │   │   └── scale_free_graph_m_*.grh    # Scale-free graphs (m=10,15,20...95)
        │   ├── graph.grh                       # Large graph file (43472 lines)
        │   ├── bn-fly-drosophila_medulla_1.grh # Drosophila brain network (1783 lines)
        │   └── cs3                             # Graph file
        │
        ├── Experimental Results/               # Output and analysis files
        │   ├── csv_output_multi_nodes_multi_desnity_*.csv # Performance benchmarks
        │   ├── output_*_log_*.txt              # Experiment logs
        │   ├── nohup.out                       # Background execution log
        │   ├── generate_synth_graphs.ipynb     # Synthetic graph generation notebook
        │   └── parse output log.ipynb          # Output analysis notebook
        │
        ├── Execution Scripts/                  # Automation and utilities
        │   ├── run.sh                          # Single algorithm execution script
        │   ├── run_both.sh                     # Comparative execution script
        │   ├── transgrh.pl                     # Graph format conversion utility
        │   └── Graph_check.py                  # Graph preprocessing utility
        │
        ├── Documentation/                      # Research papers and documentation
        │   ├── readme.md                       # Usage instructions
        │   ├── 08pce_journal.pdf              # PCE journal paper
        │   └── fpce.pdf                        # FPCE research paper
        │
        └── Build Artifacts/                    # Compiled executables and temporary files
            ├── dense-pce                       # Dense-PCE executable
            ├── dense-pce-mod                   # Modified Dense-PCE with Turan filtering
            ├── pce                             # PCE executable
            └── BBkC                            # EBBkC k-clique enumeration executable
```

## Algorithm Implementations

### 1. Dense-PCE (Original Implementation)
- **File**: `dense-pce.cpp`
- **Language**: C++
- **Purpose**: Efficient pseudo-clique enumeration with density-based optimization
- **Key Features**:
  - Uses degree-based tracking system
  - Implements sophisticated branching strategy (Type 1, 2, 3 children)
  - Supports configurable theta threshold and size constraints
  - Supports both all pseudo-cliques and maximal-only enumeration (`--maximal` flag)
  - Direct graph enumeration without external dependencies

### 2. Dense-PCE Modified (Turan-Enhanced Implementation)
- **File**: `dense-pce-mod.cpp`
- **Language**: C++
- **Purpose**: Enhanced Dense-PCE with Turan filtering for improved efficiency
- **Key Features**:
  - Integrates EBBkC for R-clique seed generation using Turan's theorem
  - Automatic R computation based on theta and minimum size parameters
  - Truss decomposition for graph preprocessing
  - Supports both all pseudo-cliques and maximal-only enumeration (`--maximal` flag)
  - Hybrid approach: uses k-clique enumeration for seeding + pseudo-clique expansion
  - Core decomposition and edge filtering optimizations

### 3. PCE (Original Implementation)
- **Directory**: `pce12/`
- **Language**: C
- **Purpose**: Original pseudo-clique enumeration algorithm with proven correctness
- **Key Features**:
  - Itemset-based approach with reverse search
  - Support for both complete enumeration (C flag) and maximal enumeration (M flag)
  - Queue-based processing with sophisticated data structures
  - Comprehensive graph format conversion utilities
  - Reference implementation for algorithm validation

### 4. FPCE (Fast Pseudo-Clique Enumeration)
- **Directory**: `FPCE/`
- **Language**: C
- **Purpose**: Optimized pseudo-clique enumeration with advanced pruning techniques
- **Key Features**:
  - Edge-bound optimization and early termination
  - Multiple algorithm variants (FPCE, PCE, ODES, FastQC)
  - Support for real-world and synthetic graph analysis
  - Biological network analysis capabilities (E. coli protein interactions)
  - ClusterONE integration for comparative analysis

### 5. EBBkC (Efficient k-clique Listing)
- **Directory**: `EBBkC/`
- **Language**: C++
- **Purpose**: High-performance k-clique enumeration with GPU acceleration
- **Key Features**:
  - Edge-oriented branching with truss decomposition
  - CUDA support for GPU acceleration
  - Early termination strategies and parallel processing
  - Comprehensive graph preprocessing and optimization
  - Integration with Dense-PCE-mod for R-clique seed generation

## Build Instructions

### Dense-PCE (Original)
```bash
cd Dense_PCE/Dense-PCE-main
g++ -std=c++17 -O3 dense-pce.cpp -o dense-pce
```

### Dense-PCE Modified (with Turan Filtering)
```bash
cd Dense_PCE/Dense-PCE-main
# First build EBBkC dependency
cd EBBkC/src && mkdir -p build && cd build
cmake .. && make
cp BBkC ../../../
cd ../../..
# Then build dense-pce-mod
g++ -std=c++17 -O3 dense-pce-mod.cpp -o dense-pce-mod -lm -lstdc++fs
```

### PCE (Original)
```bash
cd Dense_PCE/Dense-PCE-main/pce12
make
mv pce ..
```

### FPCE
```bash
cd Dense_PCE/Dense-PCE-main/FPCE/code/FPCE
make
```

### EBBkC (Standalone)
```bash
cd Dense_PCE/Dense-PCE-main/EBBkC/src
mkdir -p build && cd build
cmake ..
make
cp BBkC ../../
```

## Usage Examples

### Dense-PCE (Original)
```bash
# Enumerate all pseudo-cliques
./dense-pce synth-graphs-1000/scale_free_graph_m_10.grh --theta 0.8 --minimum 2 --maximum 100

# Enumerate only maximal pseudo-cliques
./dense-pce synth-graphs-1000/scale_free_graph_m_10.grh --theta 0.8 --minimum 2 --maximum 100 --maximal
```

### Dense-PCE Modified (with Turan Filtering)
```bash
# Enumerate all pseudo-cliques with Turan filtering (requires EBBkC)
./dense-pce-mod synth-graphs-1000/scale_free_graph_m_10.grh --theta 0.8 --minimum 3

# Enumerate only maximal pseudo-cliques with Turan filtering
./dense-pce-mod synth-graphs-1000/scale_free_graph_m_10.grh --theta 0.8 --minimum 3 --maximal
```

### PCE (Original)
```bash
# Enumerate all pseudo-cliques
./pce C -l 1 -u 100 bn-fly-drosophila_medulla_1.grh 0.9

# Enumerate only maximal pseudo-cliques
./pce M -l 1 -u 100 bn-fly-drosophila_medulla_1.grh 0.9
```

### FPCE
```bash
# Run FPCE with maximal enumeration
./FPCE/code/FPCE/fpce M -l 10 FPCE/real-synth-graph-analysis/real-graphs/bio-grid-human.txt 0.9
```

### EBBkC (for k-clique enumeration)
```bash
# Enumerate 4-cliques in Facebook dataset
./BBkC e EBBkC/dataset/facebook 4 2
```

## Graph Formats Supported

1. **EDGELIST**: Space-separated node pairs
2. **GRH**: Adjacency list format (0 to n-1 node IDs)
3. **Binary**: Optimized binary format (b_adj.bin, b_degree.bin)

## Key Features

### Graph Processing
- Multiple format support and conversion utilities
- Preprocessing for self-loops and duplicate edges
- Node ID normalization

### Algorithm Optimizations
- Degree-based tracking systems
- Early termination strategies
- Edge bounding techniques
- Parallel processing support

### Analysis Capabilities
- Real-world graph analysis
- Synthetic graph generation
- Performance benchmarking
- Memory usage tracking

### Experimental Framework
- Automated experiment scripts
- Result logging and analysis
- CSV output generation
- Jupyter notebook integration

## Dependencies

### System Requirements
- **Operating System**: Linux (Ubuntu 20.04 recommended)
- **Compiler**: GCC 9.4.0 or higher
- **Java**: OpenJDK 11.0.22 (for ClusterONE)
- **CMake**: Version 3.16+ (for EBBkC)
- **CUDA**: For GPU acceleration (optional)

### External Libraries
- **CUB**: CUDA utilities
- **ModernGPU**: GPU programming utilities
- **SparseHash**: Memory-efficient hash tables
- **Sparse++**: Sparse data structures
- **libpopcnt**: Population count operations

## Research Context

This repository implements and compares multiple approaches to pseudo-clique enumeration:

1. **Dense-PCE**: Novel density-based approach with sophisticated branching
2. **PCE**: Original pseudo-clique enumeration algorithm
3. **FPCE**: Fast variant with edge bounding optimizations
4. **EBBkC**: k-clique enumeration with early termination

The implementations support research in:
- Graph mining and analysis
- Dense subgraph discovery
- Network community detection
- Biological network analysis
- Social network analysis

## Current Status

- **Dense-PCE (Original)**: Fully implemented and functional with maximal enumeration support
- **Dense-PCE Modified**: Enhanced with Turan filtering, EBBkC integration, and maximal enumeration
- **PCE**: Original C implementation, stable reference with proven correctness
- **FPCE**: Complete with experimental framework and biological network analysis
- **EBBkC**: GPU-accelerated k-clique enumeration, integrated with Dense-PCE-mod
- **Documentation**: Comprehensive usage instructions and research paper references
- **Experiments**: Automated testing and benchmarking scripts with Jupyter notebook analysis

## Recent Enhancements

### Dense-PCE Modifications
- Added `--maximal` flag support for maximal pseudo-clique enumeration
- Implemented histogram-based maximality checking consistent with PCE.c
- Enhanced with Turan filtering via EBBkC integration for R-clique seeding
- Automatic R-value computation based on theta and minimum size parameters
- Fixed Type-3 child generation rules for accurate enumeration counts

### Integration Features
- Seamless EBBkC integration for efficient k-clique seed generation
- Hybrid enumeration approach: k-clique seeding + pseudo-clique expansion
- Cross-validation capabilities between different algorithm implementations
- Comprehensive graph preprocessing and format conversion utilities 