# Graph Pipeline Fixes and Usage Guide

## Overview

The `graph_pipeline.py` script has been significantly improved to handle path issues that occur when the project is moved to different devices or locations. The main improvements focus on robust WSL path handling and better error recovery.

## Key Fixes Applied

### 1. **Improved WSL Path Conversion**
- **Problem**: Original script didn't properly handle paths with parentheses and special characters
- **Fix**: Enhanced `_to_wsl_path()` method with single-quote path wrapping instead of escaping
- **Result**: Paths with spaces, parentheses, and special characters (like "Dense_PCE (1)") now work correctly

### 2. **Better Working Directory Handling**
- **Problem**: WSL commands didn't set the correct working directory
- **Fix**: New `_run_wsl_command()` method that properly sets working directory
- **Result**: Commands now execute in the correct directory context

### 3. **Enhanced Error Handling**
- **Problem**: Limited error reporting and no timeout handling
- **Fix**: Added comprehensive error handling with timeouts and detailed logging
- **Result**: Better debugging information and graceful failure handling

### 4. **Robust Command Execution**
- **Problem**: WSL commands could fail silently or hang
- **Fix**: Added timeout protection and better subprocess management
- **Result**: More reliable command execution with proper cleanup

## Usage

### Basic Usage

```bash
# Process a single graph file
python graph_pipeline.py testGraphs/gplus/gplus.grh

# Process with custom output directory
python graph_pipeline.py testGraphs/gplus/gplus.grh ./processed_graphs/

# Process all graphs in a directory
python graph_pipeline.py --batch synth-graphs-1000/

# Process batch with custom output base directory
python graph_pipeline.py --batch synth-graphs-1000/ ./batch_processed/
```

### Testing Path Handling

Before running the main pipeline, you can test the path handling:

```bash
python test_path_handling.py
```

This will verify:
- WSL availability and functionality
- Path conversion accuracy
- Tool accessibility
- Current directory setup

## What the Pipeline Does

The graph pipeline performs three main steps:

1. **`.grh → .edges`**: Converts the graph format to edge list format
2. **`.edges → .clean`**: Uses BBkC to clean and preprocess the edges
3. **`.clean → .bin`**: Converts to binary format for efficient processing

## Platform Support

### Windows (with WSL)
- Automatically detects WSL availability
- Uses WSL for BBkC and edgelist2binary tools
- Handles Windows path conversion to WSL paths
- Provides detailed error reporting for WSL issues

### Linux/macOS
- Uses native tools directly
- No WSL dependency
- Optimized for Unix-like systems

## Troubleshooting

### Common Issues and Solutions

#### 1. **WSL Path Conversion Errors**
```
Error: WSL command failed with return code 1
bash: line 1: cd: /mnt/c/Users/Lenovo/Downloads/Dense_PCE\ (1)/...: No such file or directory
```
**Solution**: The script now properly handles paths with parentheses and special characters. Run the test script to verify:
```bash
python test_path_handling.py
python test_path_fix.py
```

#### 2. **Missing Tools**
```
FileNotFoundError: BBkC executable not found
```
**Solution**: Ensure you're in the correct directory and tools are built:
```bash
# Check if you're in the right directory
ls -la BBkC EBBkC/

# Build tools if needed
bash build_integrated.sh
```

#### 3. **Permission Issues**
```
Permission denied when running WSL commands
```
**Solution**: Ensure WSL has proper permissions and is configured:
```bash
wsl --status
wsl --update
```

#### 4. **Timeout Errors**
```
WSL command timed out
```
**Solution**: The script now has a 5-minute timeout. For large graphs, you may need to increase this in the code.

### Debug Mode

For detailed debugging, you can modify the logging level in the script:

```python
# Change this line in graph_pipeline.py
logging.basicConfig(level=logging.DEBUG, ...)
```

## File Structure After Processing

After successful processing, you'll have:

```
input_directory/
├── graph.grh          # Original input file
├── graph.edges        # Edge list format
├── graph.clean        # Cleaned/preprocessed edges
├── b_adj.bin          # Binary adjacency data
└── b_degree.bin       # Binary degree data
```

## Performance Notes

- **Large Graphs**: The script now has timeout protection (5 minutes per step)
- **Memory Usage**: Processing is done in chunks to handle large files
- **Parallel Processing**: For batch processing, consider running multiple instances

## Integration with Build System

The graph pipeline works with the improved build system:

1. **Build the integrated system**: `bash build_integrated.sh`
2. **Process graphs**: `python graph_pipeline.py --batch synth-graphs-1000/`
3. **Run dense-pce**: `./build_integrated/dense-pce-mod-edge-order-integrated`

## Support

If you encounter issues:

1. Run the test script: `python test_path_handling.py`
2. Check the log file: `graph_preprocessing.log`
3. Verify WSL setup: `wsl --status`
4. Ensure you're in the correct directory with all required tools

The improved script should now handle location changes gracefully and provide better error reporting when issues occur.
