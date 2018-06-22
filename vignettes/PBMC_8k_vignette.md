---
title: "Doublet Detection on 8k PBMCs from 10X Genomics"
author: "Tom Kelly and Jonathon Shor"
date: "2018-06-22"
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
#> Downloading GitHub repo TomKellyGenetics/Rphenograph@master
#> from URL https://api.github.com/repos/TomKellyGenetics/Rphenograph/zipball/master
#> Installing Rphenograph
#> '/opt/local/R-3.5/lib/R/bin/R' --no-site-file --no-environ --no-save  \
#>   --no-restore --quiet CMD INSTALL  \
#>   '/tmp/RtmpDpAijD/devtools6e2a30cb85fa/TomKellyGenetics-Rphenograph-456c59a'  \
#>   --library='/home/tom/R/x86_64-pc-linux-gnu-library/3.5'  \
#>   --install-tests
#> 
#> '/opt/local/R-3.5/lib/R/bin/R' --no-site-file --no-environ --no-save  \
#>   --no-restore --quiet CMD INSTALL  \
#>   '/tmp/RtmpDpAijD/devtools6e2a39c19e7e/TomKellyGenetics-DoubletDetection-c8c79b5'  \
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
#> Loading required package: Matrix
#> Loading required package: Rphenograph
#> Loading required package: ggplot2
#> Loading required package: igraph
#> 
#> Attaching package: 'igraph'
#> The following objects are masked from 'package:stats':
#> 
#>     decompose, spectrum
#> The following object is masked from 'package:base':
#> 
#>     union
#> Error: package or namespace load failed for 'Rphenograph' in library.dynam(lib, package, package.lib):
#>  shared object 'Rphenograph.so' not found
#> Error: package 'Rphenograph' could not be loaded
```

## Download Data from 10X Genomics

These commands will be called from R in Linux and Mac systems. The files will need to be downloaded from the given url (into the working directory) in Windows.


```r
#download files
if(!(file.exists("pbmc8k_filtered_gene_bc_matrices.tar.gz"))){
  system("wget http://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc8k/pbmc8k_filtered_gene_bc_matrices.tar.gz")
  #extract compressed files 
  #system("tar -xvzf pbmc8k_filtered_gene_bc_matrices.tar.gz")
}
```

## Load Counts Matrix


```r
matrix_path <- 'filtered_gene_bc_matrices/GRCh38/matrix.mtx'
raw_counts <- DoubletDetection::load_mtx(matrix_path)
#> Error in library.dynam(lib, package, package.lib): shared object 'Rphenograph.so' not found
# Remove columns with all 0s
non_zero_genes <- apply(raw_counts, 1, function(x) sum(x) != 0)
#> Error in apply(raw_counts, 1, function(x) sum(x) != 0): object 'raw_counts' not found
raw_counts <- raw_counts[non_zero_genes, ]
#> Error in eval(expr, envir, enclos): object 'raw_counts' not found
```

## Run Doublet Detection

Right now, phenograph is a bit talkative, so we suppress the output to avoid lots of text:


```
#> Error in library.dynam(lib, package, package.lib): shared object 'Rphenograph.so' not found
#> Error in eval(expr, envir, enclos): object 'clf' not found
```


```r
print(paste("Time elapsed", end-start,"seconds,", (end-start) / clf.n_iters, clf.n_iters), "sec/iteration, for", n_inters, "iterations")
#> Error in paste("Time elapsed", end - start, "seconds,", (end - start)/clf.n_iters, : object 'clf.n_iters' not found
```

## Visualize Results

### Convergence of Doublet Calls


```r
DoubletDetection::convergence(clf, save='convergence_test.pdf', show=TRUE)
#> Error in library.dynam(lib, package, package.lib): shared object 'Rphenograph.so' not found
```


### Doublets on tSNE map


```r
tsne_coords <- DoubletDetection::tsne_plot(raw_counts, doublets, save='tsne_test.pdf', show=True)
#> Error in library.dynam(lib, package, package.lib): shared object 'Rphenograph.so' not found
```




