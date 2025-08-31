#!/usr/bin/env python3
"""
Test script for path handling in the graph pipeline.
This script tests the WSL path conversion and command execution.
"""

import sys
import os
from pathlib import Path
import subprocess
import platform

def test_wsl_path_conversion():
    """Test the WSL path conversion function."""
    print("=== Testing WSL Path Conversion ===")
    
    # Test cases with different path types
    test_paths = [
        Path("C:/Users/Test User/Documents/file.txt"),
        Path("C:/Users/Lenovo/Downloads/Dense_PCE (1)/Dense_PCE/CSE491-/Dense_PCE/Dense-PCE-main"),
        Path("D:/OneDrive - northsouth.edu/Desktop/Dense_PCE/CSE491-/Dense_PCE/Dense-PCE-main"),
        Path("E:/Dense_PCE/CSE491-/Dense_PCE/Dense-PCE-main"),
        Path("/home/user/test"),
        Path("relative/path/file.txt")
    ]
    
    for test_path in test_paths:
        print(f"Original: {test_path}")
        # Simulate the fixed _to_wsl_path function
        p = str(test_path.resolve())
        if len(p) > 1 and p[1] == ':' and p[0].isalpha():
            drive = p[0].lower()
            rest = p[2:].replace('\\', '/')
            # Don't escape here - we'll use single quotes around the path
            wsl_path = f"/mnt/{drive}{rest}"
        else:
            wsl_path = p.replace('\\', '/')
        
        print(f"WSL Path: {wsl_path}")
        print()

def test_wsl_availability():
    """Test if WSL is available and working."""
    print("=== Testing WSL Availability ===")
    
    is_windows = platform.system() == "Windows"
    print(f"Platform: {platform.system()}")
    print(f"Is Windows: {is_windows}")
    
    if is_windows:
        # Check if wsl command is available
        try:
            result = subprocess.run(["wsl", "--version"], capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                print("✓ WSL is available")
                print(f"WSL version: {result.stdout.strip()}")
            else:
                print("✗ WSL command failed")
        except FileNotFoundError:
            print("✗ WSL command not found")
        except subprocess.TimeoutExpired:
            print("✗ WSL command timed out")
        except Exception as e:
            print(f"✗ Error checking WSL: {e}")
    else:
        print("Not on Windows, WSL not applicable")

def test_current_directory():
    """Test current working directory and path resolution."""
    print("=== Testing Current Directory ===")
    
    cwd = Path.cwd()
    print(f"Current working directory: {cwd}")
    print(f"Absolute path: {cwd.absolute()}")
    print(f"Resolved path: {cwd.resolve()}")
    
    # Check if we're in the expected project directory
    if (cwd / "dense-pce-mod-edge-order.cpp").exists():
        print("✓ Found dense-pce-mod-edge-order.cpp - likely in correct directory")
    else:
        print("✗ dense-pce-mod-edge-order.cpp not found - may be in wrong directory")
    
    # Check for key directories
    key_dirs = ["EBBkC", "BBkC", "testGraphs"]
    for dir_name in key_dirs:
        if (cwd / dir_name).exists():
            print(f"✓ Found {dir_name}/")
        else:
            print(f"✗ {dir_name}/ not found")

def test_tool_paths():
    """Test if required tools are accessible."""
    print("=== Testing Tool Paths ===")
    
    cwd = Path.cwd()
    
    # Check BBkC
    bbkc_path = cwd / "BBkC"
    if bbkc_path.exists():
        print(f"✓ BBkC found at: {bbkc_path}")
        if bbkc_path.is_file():
            print("  - BBkC is a file (executable)")
        else:
            print("  - BBkC is a directory")
    else:
        print(f"✗ BBkC not found at: {bbkc_path}")
    
    # Check edgelist2binary
    edgelist2binary_path = cwd / "EBBkC" / "Cohesive_subgraph_book" / "datasets" / "edgelist2binary"
    if edgelist2binary_path.exists():
        print(f"✓ edgelist2binary found at: {edgelist2binary_path}")
        if edgelist2binary_path.is_file():
            print("  - edgelist2binary is a file (executable)")
        else:
            print("  - edgelist2binary is a directory")
    else:
        print(f"✗ edgelist2binary not found at: {edgelist2binary_path}")

def test_wsl_command_execution():
    """Test WSL command execution with a simple command."""
    print("=== Testing WSL Command Execution ===")
    
    if platform.system() != "Windows":
        print("Not on Windows, skipping WSL test")
        return
    
    try:
        # Test a simple WSL command
        result = subprocess.run(
            ["wsl", "bash", "-c", "pwd && echo 'WSL is working'"],
            capture_output=True,
            text=True,
            timeout=10
        )
        
        if result.returncode == 0:
            print("✓ WSL command execution successful")
            print(f"Output: {result.stdout.strip()}")
        else:
            print("✗ WSL command failed")
            print(f"Error: {result.stderr}")
            
    except Exception as e:
        print(f"✗ Error testing WSL command execution: {e}")

def main():
    """Run all tests."""
    print("Path Handling Test Suite")
    print("=" * 50)
    
    test_current_directory()
    print()
    
    test_tool_paths()
    print()
    
    test_wsl_availability()
    print()
    
    test_wsl_path_conversion()
    print()
    
    test_wsl_command_execution()
    print()
    
    print("Test suite completed!")

if __name__ == "__main__":
    main()
