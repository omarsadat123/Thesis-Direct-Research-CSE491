# FPCE: Optimized Exact Pseudo-Clique Enumeration

This repository collects the research artifacts for the CSE491 direct research project on **FPCE**, with supporting material for **DensePCE** and related pseudo-clique enumeration experiments. The main thesis report is the canonical reference for methodology, implementation details, and results.

## Table of Contents

- [Summary](#summary)
- [Repository Structure](#repository-structure)
- [Thesis Reference](#thesis-reference)
- [Quick Start](#quick-start)
- [Citation](#citation)

## Summary

The project combines:

- notebook-based experimentation and analysis
- C and C++ implementations for enumeration and bound-based experiments
- datasets used for testing and benchmarking
- the final thesis report in PDF format

## Repository Structure

| Path | Purpose |
| --- | --- |
| `Dense_PCE/` | Top-level area for DensePCE work, notebooks, and experiment artifacts. |
| `Dense_PCE/Dense-PCE-main/` | Main DensePCE implementation and research workspace. |
| `Dense_PCE/Dense-PCE-main/pce12/` | PCE-related code and experiment-specific files. |
| `OrderBound/` | Order-bound implementation and supporting C/C++ research code. |
| `dataset/` | Graph benchmark datasets in `.grh` format. |
| `Graph/` | Additional graph inputs, examples, or supporting materials. |
| `Thesis report direct research.pdf` | Primary thesis document and canonical project reference. |

## Thesis Reference

The thesis report is the best starting point for understanding the project goals, experimental setup, and conclusions:

- [`Thesis report direct research.pdf`](./Thesis%20report%20direct%20research.pdf)

## Quick Start

A minimal workflow consistent with the repository contents:

1. Read the thesis PDF to understand the research context.
2. Inspect the notebooks in `Dense_PCE/` for experimentation and analysis.
3. Review the implementation code in `OrderBound/` and `Dense_PCE/Dense-PCE-main/`.
4. Explore the datasets in `dataset/` and `Graph/` for benchmarking or reproduction.

## Citation

If you use this repository in academic work, please cite the thesis appropriately.
