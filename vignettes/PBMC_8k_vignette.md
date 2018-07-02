---
title: "Doublet Detection on 8k PBMCs from 10X Genomics"
author: "Tom Kelly and Jonathon Shor"
date: "2018-06-29"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Vignette Title}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



## Install packages

### Install package:


```r
if(!require(devtools)){
  install.packages("devtools") # If not already installed
}
#> Loading required package: devtools
devtools::install_github("TomKellyGenetics/DoubletDetection", ref = "r-implementation")
#> Downloading GitHub repo TomKellyGenetics/DoubletDetection@r-implementation
#> from URL https://api.github.com/repos/TomKellyGenetics/DoubletDetection/zipball/r-implementation
#> Installing DoubletDetection
#> '/opt/local/R-3.5/lib/R/bin/R' --no-site-file --no-environ --no-save  \
#>   --no-restore --quiet CMD INSTALL  \
#>   '/tmp/RtmpuYkjNR/devtools30368d78cbe/TomKellyGenetics-DoubletDetection-05ab87f'  \
#>   --library='/home/tom/R/x86_64-pc-linux-gnu-library/3.5'  \
#>   --install-tests
#> 
#> Installation failed: Command failed (1)
```

### Install dependancies:







```r
#devtools::install_github("JinmiaoChenLab/Rphenograph")
devtools::install_github("TomKellyGenetics/Rphenograph")
install.packages("hdf5r")
install.packages("Matrix")
```

### Install hdf5 (if needed)

OS X (using Homebrew)	`brew install hdf5`

Debian-based systems (including Ubuntu)	`sudo apt-get install libhdf5-dev`

Systems supporting yum and RPMs	`sudo yum install hdf5-devel`

HDF5 1.8.14 has been pre-compiled for Windows and is available at https://github.com/mannau/h5-libwin 

### load package (and dependancies)


```r
library("DoubletDetection")
#> Error in library("DoubletDetection"): there is no package called 'DoubletDetection'
```

## Download Data from 10X Genomics

These commands will be called from R in Linux and Mac systems. The files will need to be downloaded from the given url (into the working directory) in Windows.


```r
#download files
if(!(file.exists("pbmc8k_filtered_gene_bc_matrices.tar.gz"))){
  system("wget http://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc8k/pbmc8k_filtered_gene_bc_matrices.tar.gz")
  #extract compressed files 
  system("tar -xvzf pbmc8k_filtered_gene_bc_matrices.tar.gz")
}
```

## Load Counts Matrix


```r
matrix_path <- 'filtered_gene_bc_matrices/GRCh38/matrix.mtx'
raw_counts <- DoubletDetection::load_mtx(matrix_path)
#> Error in loadNamespace(name): there is no package called 'DoubletDetection'
# Remove columns with all 0s
non_zero_genes <- apply(raw_counts, 1, function(x) sum(x) != 0)
#> Error in apply(raw_counts, 1, function(x) sum(x) != 0): object 'raw_counts' not found
raw_counts <- raw_counts[non_zero_genes, ]
#> Error in eval(expr, envir, enclos): object 'raw_counts' not found
```

## Run Doublet Detection

Right now, phenograph is a bit talkative, so we suppress the output to avoid lots of text:












