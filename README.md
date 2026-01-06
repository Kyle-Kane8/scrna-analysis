# Single-Cell RNA-seq Analysis: Macrophage & Immune Profiling

This repository contains an end-to-end single-cell RNA-seq workflow using Scanpy.

It includes:
- Synthetic 10X-style count matrix (for demonstration only)
- QC, normalization, dimensionality reduction
- Louvain/Leiden clustering
- Marker identification and visualization

## Synthetic Data

A small synthetic 10X-format dataset is provided under `data/`:
- `matrix.mtx`
- `genes.tsv`
- `barcodes.tsv`

This allows the pipeline to run without downloading external data, but the
values are randomly generated and not biologically meaningful.
