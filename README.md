# scClustering-Benchmark

# From general metrics to local practice: a comprehensive benchmark of clustering methods for single-cell transcriptomics

This repository contains the analysis code associated with the manuscript:

**From general metrics to local practice: a comprehensive benchmark of clustering methods for single-cell transcriptomics**

Authors: Ziying Huang, Zedong Lin, Anqi Zou, Jiaqi Xu, Yue Jiang, Elaine Huang, Yue Wang, Chenxi Xiao, Yunshu Dong, Jou-Hsuan Lee, Hongjuan Niu, Ping Zhou, Rong Li, Chao Zhang, and Heng Pan.

## Overview

Clustering is a central step in single-cell RNA sequencing (scRNA-seq) analysis and directly affects downstream interpretation, including cell-type annotation, differential expression, trajectory analysis, and biological discovery. This study benchmarks representative clustering methods across diverse real-world scRNA-seq datasets and evaluates their performance from three complementary perspectives:

1. **General clustering performance** across public benchmark datasets.
2. **Local performance in biologically challenging scenarios**, with a focus on endometrial datasets involving subtle transcriptional variation and residual low-quality cells.
3. **Practical usability**, including installation, documentation, workflow compatibility, parameter transparency, and reproducibility support.

The benchmark includes **16 clustering methods**, **36 publicly available scRNA-seq datasets**, and **5 clustering evaluation metrics**.

## Benchmark design

The analysis was designed to evaluate clustering methods under both broad and context-specific conditions.

### Clustering methods

The benchmark includes both conventional and deep learning-based methods.

**Non-deep learning methods**

- CIDR
- SC3
- RaceID
- VarID
- SHARP
- Spectrum
- Seurat
- Monocle3
- TSCAN
- Scanpy

**Deep learning-based methods**

- scGAC
- scDeepCluster
- scDCC
- scTAG
- DeepScena

For scDeepCluster, both manually specified cluster number and automatic cluster number selection were evaluated.

Detailed package versions and software sources are summarized in `Supplementary Table 1`.

### Datasets

The study uses publicly available scRNA-seq datasets covering different species, tissues, protocols, cell numbers, and annotation schemes. These datasets include embryonic, immune, neural, cardiac, intestinal, reproductive, and developmental systems, with a specific emphasis on human endometrial datasets for biologically challenging local evaluation.

Detailed dataset information, including accession number, organism, tissue, protocol, cell number, feature number, cluster number, and label type, is summarized in `Supplementary Table 2`.

### Evaluation metrics

Clustering performance was evaluated using five complementary metrics:

- Adjusted Rand Index (ARI)
- Fowlkes-Mallows Index (FMI)
- Adjusted Mutual Information (AMI)
- Normalized Mutual Information (NMI)
- Homogeneity

Detailed descriptions of these metrics are provided in `Supplementary Table 3`.

## Data availability

This repository does not contain large raw sequencing files or complete processed single-cell objects. Raw and processed data should be downloaded from the public repositories listed in the manuscript and Supplementary Table 2.

Datasets used in this benchmark include, but are not limited to:

- E-MTAB-3321
- GSE65525
- GSE189120
- GSE171993
- GSE60361
- PBMC CITE-seq dataset from the scvi-tools benchmark resource
- fig11981034
- GSE137537
- GSE109816
- E-MTAB-6701
- GSE125970
- GSE109774
- GSE201333
- GSE130664
- GSE111976
- E-MTAB-10287
- GSE119945

After download, each dataset should be organized into the following standardized format before running the benchmark:

```text
data/
└── processed_10x/
    └── <dataset_name>/
        ├── matrix.mtx.gz
        ├── features.tsv.gz
        ├── barcodes.tsv.gz
        └── labels.csv
```

The `labels.csv` file should contain the ground-truth cell-type annotation used for benchmarking.

## Software requirements

The analyses were performed using R and Python. Major software and packages include:

### R

- R
- Seurat
- monocle3
- SC3
- RaceID
- SHARP
- Spectrum
- TSCAN
- Matrix
- SingleCellExperiment
- dplyr
- tidyr
- ggplot2
- patchwork
- pheatmap
- ComplexHeatmap
- mclust
- aricode

### Python

- Python
- scanpy
- anndata
- numpy
- pandas
- scipy
- scikit-learn
- igraph
- leidenalg
- matplotlib
- seaborn
- torch, if required by deep learning-based methods

Because some deep learning-based methods require method-specific environments, we recommend creating separate conda environments for each deep learning tool when necessary.

Full package versions should be recorded in:

```text
docs/session_info.txt
docs/software_versions.txt
```

## License

This repository is released under the MIT License. See `LICENSE` for details.

## Contact

For questions about the code or benchmark framework, please contact the corresponding authors listed in the manuscript.

## Notes

This repository is intended to provide the custom scripts required to reproduce the analyses in the manuscript. It does not redistribute large public sequencing datasets. Users should download raw or processed data from the original public repositories and organize them according to the instructions above before running the analysis.
