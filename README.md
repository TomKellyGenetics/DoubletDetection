#Notice

**Do not attempt to run code from this package**

**This repository is still in progress**

The R implementation of functions provided in the Python version is ongoing. This package is not functional yet.

Contributions to the R implementation are welcome.

# DoubletDetection

DoubletDetection is an R implementation of a package to detect doublets (technical errors) in single-cell RNA-seq count matrices.

## Installation

To install DoubletDetection in R:

```
if(!require(devtools)){
  install.packages("devtools") # If not already installed
}
devtools::install_github("TomKellyGenetics/DoubletDetection", ref = "r-implementation")
```

## Running

To run basic doublet classification:

```
library("DoubletDetection")
clf <- BoostClassifier$new()
# raw_counts is a cells by genes count matrix
labels = clf$fit(raw_counts)$predict()
```

`raw_counts` is a scRNA-seq count matrix (cells by genes) or data.frame. `labels` is a binary numerical vector with the value 1 representing a 
detected doublet.

# Usage

## R Version

These functions and methods (for the Reference Class) have been documented and can be accessed in the R help system. A [vignette][https://rawgit.com/TomKellyGenetics/DoubletDetection/r-implementation/vignettes/PBMC_8k_vignette.html] will be prepared using Rmarkdown in due course.

## Python Version

For the Python implementation, see the original repository: https://github.com/JonathanShor/DoubletDetection

See their [jupyter notebook](https://nbviewer.jupyter.org/github/JonathanShor/DoubletDetection/blob/master/docs/PBMC_8k_vignette.ipynb) for an example on 8k PBMCs from 10x.

## Obtaining data
Data can be downloaded from the [10x website](https://support.10xgenomics.com/single-cell/datasets).


## Citations

bioRxiv submission is in the works.

This project is licensed under the terms of the MIT license.
