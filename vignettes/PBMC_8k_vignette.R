if(!require(devtools)){
  install.packages("devtools") # If not already installed
}
devtools::install_github("TomKellyGenetics/DoubletDetection", ref = "r-implementation")

### Install dependancies:

options(repos = "https://cran.stat.auckland.ac.nz/")
#devtools::install_github("JinmiaoChenLab/Rphenograph")
devtools::install_github("TomKellyGenetics/Rphenograph")
#install.packages("hdf5r")
#install.packages("Matrix")

library("DoubletDetection")

#download files
if(!(file.exists("pbmc8k_filtered_gene_bc_matrices.tar.gz"))){
  system("wget http://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc8k/pbmc8k_filtered_gene_bc_matrices.tar.gz")
  #extract compressed files 
  #system("tar -xvzf pbmc8k_filtered_gene_bc_matrices.tar.gz")
}

matrix_path <- 'filtered_gene_bc_matrices/GRCh38/matrix.mtx'
raw_counts <- DoubletDetection::load_mtx(matrix_path)
# Remove columns with all 0s
non_zero_genes <- apply(raw_counts, 1, function(x) sum(x) != 0)
raw_counts <- raw_counts[non_zero_genes, ]

#initialise BoostClassifier Reference Class (with analysis parameters)
clf <- DoubletDetection::BoostClassifier$new(n_iters=50)
 
#run doublet prediction
start <- Sys.time()
doublets <- clf$fit(raw_counts)$predict()
end <- Sys.time()

save.image("vignette_test_source.RData")

print(paste("Time elapsed", end-start,"seconds,", (end-start) / clf.n_iters, clf.n_iters), "sec/iteration, for", n_inters, "iterations")

DoubletDetection::convergence(clf, save='convergence_test.pdf', show=TRUE)

tsne_coords <- DoubletDetection:tsne_plot(raw_counts, doublets, save='tsne_test.pdf', show=True)


