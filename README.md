# Notice

**This repository is still in development**

The testing of R implementation of functions provided in the Python version is ongoing. This package is not fully documented yet. Please share any problems you encounter in the documentation or functionality of the code as an issue on GitHub. Thank you for your patience.

Contributions to the R implementation via pull request are welcome.

# DoubletDetection

## Version 2.3.0

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/DoubletDetection)](https://cran.r-project.org/package=DoubletDetection)
[![Travis Build Status](https://travis-ci.org/TomKellyGenetics/DoubletDetection.svg?branch=r-implementation)](https://travis-ci.org/TomKellyGenetics/DoubletDetection)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/TomKellyGenetics/DoubletDetection?branch=r-implementation&svg=true)](https://ci.appveyor.com/project/TomKellyGenetics/DoubletDetection)
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![codecov](https://codecov.io/gh/TomKellyGenetics/DoubletDetection/branch/r-implementation/graph/badge.svg)](https://codecov.io/gh/TomKellyGenetics/DoubletDetection)


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
# raw_counts is a genes by cells count matrix
labels = clf$fit(raw_counts)$predict()
```

`raw_counts` is a scRNA-seq count matrix or data.frame (genes x cells).
Note that the dimensions of input matrix differs from the Python version.

 `labels` is a binary numerical vector with the value `1` representing a 
detected doublet, `0` a singlet, and `NA` an ambiguous cell.

The classifier works best when there are several cell types present in the data. Furthermore, it should be applied individually to each run in an aggregated count matrix.

# Usage

## R Version

These functions and methods (for the Reference Class) have been documented and can be accessed in the R help system. A [vignette](https://rawgit.com/TomKellyGenetics/DoubletDetection/r-implementation/vignettes/PBMC_8k_vignette.html) will be prepared using Rmarkdown in due course.

## Python Version

For the Python implementation, see the original repository: https://github.com/JonathanShor/DoubletDetection

See their [jupyter notebook](https://nbviewer.jupyter.org/github/JonathanShor/DoubletDetection/blob/master/docs/PBMC_8k_vignette.ipynb) for an example on 8k PBMCs from 10x.

## Obtaining data
Data can be downloaded from the [10x website](https://support.10xgenomics.com/single-cell/datasets).


# Citations

## R version

Please cite the R implementation as an R package using `citation(DoubletDetection)`.

>Adam J. Gayoso, Jonathan D. Shor, Ambrose J. Carr, and S. Thomas Kelly (2018). DoubletDetection: a package to detect
doublets (technical errors) in single-cell RNA-seq count matrices. R package version 2.3.0
https://github.com/TomKellyGenetics/DoubletDetection"

## Python version

Please acknowledge the original contributors when using the R implementation.

bioRxiv submission is in progress. Please refer to the [Python Repository]() https://github.com/JonathanShor/DoubletDetection) for more details.

This project is licensed under the terms of the MIT license (in accordance with the license of the original repository).
