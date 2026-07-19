# FPCE*: An Optimized Exact Pseudo-Clique Enumerator

This repository contains the implementation artifacts and experimental materials for the directed research thesis **“FPCE*: An Optimized Exact Pseudo-Clique Enumerator”** (CSE491/ECE, North South University).

The work focuses on exact \((\ell, \theta)\)-pseudo-clique enumeration on large graphs, with implementation and optimization choices grounded in the thesis report.

## Table of Contents

- [Project Summary](#project-summary)
- [Repository Structure](#repository-structure)
- [Thesis Report (Primary Reference)](#thesis-report-primary-reference)
- [Quick Start](#quick-start)
- [Citation](#citation)

## Project Summary

According to the thesis report, the project builds an optimized exact pseudo-clique enumerator around the DensePCE/FPCE line of work, including:

- reverse-search based enumeration
- pruning strategies such as order-bound and core-aware bounds
- optimized graph/data handling for large real-world datasets
- benchmarking on multiple graph datasets and synthetic graphs

For full methodology, proofs, and result analysis, use the thesis PDF linked below.

## Repository Structure

| Path | Purpose |
| --- | --- |
| `Dense_PCE/` | Main Dense-PCE research workspace (notebooks, scripts, C/C++ code, and PCE sources). |
| `Dense_PCE/Dense-PCE-main/` | Core experimentation area (`dense-pce.cpp`, notebooks, runner scripts, generated outputs). |
| `Dense_PCE/Dense-PCE-main/pce12/` | Upstream/legacy PCE implementation and build files. |
| `OrderBound/` | Order-bound related implementation variant (`dense-pce.cpp`). |
| `dataset/` | Real graph benchmark datasets in `.grh` format used in experiments. |
| `Graph/` | Additional graph data bundle (`graph_dataset.zip`). |
| `Thesis report direct research.pdf` | Full thesis document with problem definition, methodology, experiments, and conclusions. |

## Thesis Report (Primary Reference)

Read the thesis first for complete context:

- [`Thesis report direct research.pdf`](./Thesis%20report%20direct%20research.pdf)

This is the authoritative source for:

- formal definitions (clique, edge density, \((\ell,\theta)\)-pseudo-clique)
- algorithmic design and optimization rationale
- experimental setup, ablation details, and performance comparisons

## Quick Start

### 1) Build the original PCE binary

```bash
cd Dense_PCE/Dense-PCE-main/pce12
make
```

### 2) Build Dense-PCE implementation

```bash
cd Dense_PCE/Dense-PCE-main
g++ -O3 dense-pce.cpp -o dense-pce
```

### 3) Explore inputs and experiments

- Datasets: `dataset/*.grh`
- Synthetic graph generation: `Dense_PCE/Dense-PCE-main/generate_synth_graphs.ipynb`
- Existing run helpers: `Dense_PCE/Dense-PCE-main/run.sh`, `run_both.sh`

## Citation

If you use this repository in academic work, please cite the thesis report:

- FPCE*: An Optimized Exact Pseudo-Clique Enumerator (Directed Research, North South University, Spring 2026).
