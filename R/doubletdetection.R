##' @name doublet_detection
##' @rdname normalize_counts
##'
##' @title Doublet detection in single-cell RNA-seq data
##' 
##' @description Functions required to compute identification of doublets in single-cell RNA-Seq experiments.
##' 
##' @param raw_counts a matrix or data.frame containing raw UMI counts such as a gene-barcode matrix  (genes or transcripts by cell barcode).
##' @param pseudocount a numeric to add prior to log-transform (to avoid non-zero values).
##' 
##' @keywords scRNA quality filter matrix normalise doublets single-cell
##' @export
normalize_counts <- function(raw_counts, pseudocount=0.1){
  # Normalize count array.
  # 
  # Args:
  # raw_counts (ndarray): count data
  # pseudocount (numeric, optional): Count to add prior to log transform.
  # 
  # Returns:
  # ndarray: Normalized data.
  # Sum across cells
  cell_sums <- apply(raw_counts, 2, sum)
  
  # Mutiply by median and divide each cell by cell sum
  median <- median(cell_sums)
  normed <- apply(raw_counts, 2, function(x) x*median/sum(x))
  
  #log-transform
  normed <- log(normed + pseudocount)
  
  return(normed)
}


##' @rdname load_mtx
##' 
##' @title Load count matrix in mtx format
##' 
##' @param file character: path to mtx file
##' @import Matrix
##' 
##' @export
load_mtx <- function(file){
  #     Load count matrix in mtx format
  # 
  #     Args:
  #         file (str): Path to mtx file
  # 
  #     Returns:
  #         ndarray: Raw count matrix.
  if(file.exists(file)){
    group <- normalizePath(paste(dirname(file)))
  } else{
    warning("That genome does not exist in this file")
    return(NULL)
  }
  gene_names <- read.table(paste(group, "genes.tsv", sep = "/"))[,1]
  barcodes <- readLines(paste(group, "barcodes.tsv", sep = "/"))
  raw_counts <- readMM(paste(group, "matrix.mtx", sep = "/"))
  rownames(raw_counts) <- gene_names
  colnames(raw_counts) <- barcodes
  raw_counts <- as.matrix(raw_counts)
  return(raw_counts)
}

##' @rdname load_10x_h5
##' 
##' @title  Load count matrix in 10x H5 format
##' 
##' @description Adapted from: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices
##' 
##' @param file character: path to H5 file
##' @param genome character: path to top level H5 group
##' @import hdf5r Matrix
##' 
##' @export
load_10x_h5 <- function(file, genome = NULL, barcode_filtered = TRUE){
  #    
  #        Adapted from:
  #        https://support.10xgenomics.com/single-cell-gene-expression/software/
  #        pipelines/latest/advanced/h5_matrices
  # 
  #     Args:
  #         file (str): Path to H5 file
  #         genome (str): genome, top level h5 group
  # 
  #     Returns:
  #         ndarray: Raw count matrix.
  #         ndarray: Barcodes
  #         ndarray: Gene names
  if (!file.exists(normalizePath((file)))) 
    stop("Could not find the path of file: '", file, 
         "'. Please double-check if the directory exists.\n")
  if (!dir.exists(normalizePath(paste(dirname(file))))) 
    stop("Could not find the pipestance output directory: '", 
         normalizePath(paste(dirname(file))), "'. Please double-check if the directory exists.\n")
  h5_path <- file.path(normalizePath(paste(dirname(file))), ifelse(barcode_filtered, 
                                                                   "filtered_gene_bc_matrices_h5.h5", "raw_gene_bc_matrices_h5.h5"))
  if (!file.exists(h5_path)) {
    stop(sprintf("Could not find matrix H5: %s\n", h5_path))
  }
  file.h5 <- H5File$new(h5_path, mode = "r")
  #names(file.h5)
  #file.h5$ls(recursive=TRUE)
  
  # gene_ids <- file.h5[[paste0(genome, "/genes")]][]
  #gene_names <- eval(parse(text = paste0("file.h5[[\"", genome, "\"]][[\"gene_names\"]]")))[]
  #barcodes <- eval(parse(text = paste0("file.h5[[\"", genome, "\"]]")))[["barcodes"]]
  gene_names <- file.h5[[paste0(genome, "/gene_names")]][]
  barcodes <- file.h5[[paste0(genome, "/barcodes")]][]
  data <- file.h5[[paste0(genome, "/data")]][]
  indices <- file.h5[[paste0(genome, "/indices")]][]
  indptr <- file.h5[[paste0(genome, "/indptr")]][]
  shape <- file.h5[[paste0(genome, "/shape")]][]
  #h5attr_names(file.h5[[paste0(genome, "/shape")]])
  matrix <- sparseMatrix(i = indices+1, p = indptr, x = data, dims = shape)
  dense_matrix <- matrix(matrix, nrow = shape[1], ncol = shape[2])
  
  rownames(dense_matrix) <- gene_names
  colnames(dense_matrix) <- barcodes
  
  return(dense_matrix)
}

##' @rdname BoostClassifier
##' @title Classifier for doublets in single-cell RNA-seq data
##' 
##' @import Matrix stats Rphenograph
##' 
##' @export
##' @field boost_rate (numeric, optional): Proportion of cell population size to produce as synthetic doublets.
##' @field n_components (integer optional): Number of principal components used for clustering.
##' @field n_top_var_genes (integer optional): Number of highest variance genes to use; other genes discarded. Will use all genes when zero.
##' @field new_lib_as: (([integer integer]) -> integer optional): Method to use in choosing library size for synthetic doublets. Defaults to NULL which makes synthetic doublets the exact addition of its parents; alternative is new_lib_as = max.
##' @field replace (logical, optional): If FALSE, a cell will be selected as a synthetic doublet's parent no more than once.
##' @field phenograph_parameters (list, optional): Parameter dict to pass directly to Phenograph. Note that we change the Phenograph 'prune' default to TRUE; you must specifically include 'prune': FALSE here to change this.
##' @field n_iters (integer optional): Number of fit operations from which to collect p-values. Defualt value is 25. normalizer ((matrix) -> matrix): Method to normalize raw_counts. Defaults to normalize_counts, included in this package. Note: To use normalize_counts with its pseudocount parameter changed from the default 0.1 value to some positive numeric `new_var`, use: normalizer=lambda counts: doubletdetection.normalize_counts(counts, pseudocount=new_var)
##' @field normalizer ((matrix) -> matrix): Method to normalize raw_counts. Defaults to normalize_counts, included in this package. Note: To use normalize_counts with its pseudocount parameter changed from the default 0.1 value to some positive numeric `new_var`, use: normalizer=lambda counts: doubletdetection.normalize_counts(counts, pseudocount=new_var)
##' 
##' @param p_thresh (numeric, optional): hypergeometric test p-value threshold that determines per iteration doublet calls
##' @param voter_thresh (numeric, optional): fraction of iterations a cell must be called a doublet
##' @param cell1,cell2 (vector, numeric): Gene count vectors.

BoostClassifier <- setRefClass(
  "BoostClassifier",
  fields = list(
    boost_rate = "numeric", #ANY for no type specified (all passed on)
    n_components = "integer",
    n_top_var_genes = "integer",
    new_lib_as = "ANY",
    replace = "logical",
    phenograph_parameters = "list",
    n_iters = "integer",
    normalizer = "function"#,
    #raw_counts = "ANY",
    #p_thres = "numeric",
    #voter_thres = "numeric",
    #cell1 = "numeric",
    #cell2 = "numeric"
  ),
  methods = list(
    initialize = function(boost_rate = 0.25,
                          n_components = 30L,
                          n_top_var_genes = 10000L,
                          new_lib_as = NULL,
                          replace = FALSE, 
                          phenograph_parameters = list(prune = TRUE),
                          n_iters = 25L,
                          normalizer = normalize_counts){
      "This method is called when you create an instance of the class."
      #     Parameters:
      #         boost_rate (numeric, optional): Proportion of cell population size to
      #             produce as synthetic doublets.
      #         n_components (integer optional): Number of principal components used for
      #             clustering.
      #         n_top_var_genes (integer optional): Number of highest variance genes to
      #             use; other genes discarded. Will use all genes when zero.
      #         new_lib_as: (([integer integer]) -> integer optional): Method to use in choosing
      #             library size for synthetic doublets. Defaults to NULL which makes
      #             synthetic doublets the exact addition of its parents; alternative
      #             is new_lib_as = max.
      #         replace (logical, optional): If FALSE, a cell will be selected as a
      #             synthetic doublet's parent no more than once.
      #         phenograph_parameters (dict, optional): Parameter dict to pass directly
      #             to Phenograph. Note that we change the Phenograph 'prune' default to
      #             TRUE; you must specifically include 'prune': FALSE here to change
      #             this.
      #         n_iters (integer optional): Number of fit operations from which to collect
      #             p-values. Defualt value is 25.
      #         normalizer ((ndarray) -> ndarray): Method to normalize raw_counts.
      #             Defaults to normalize_counts, included in this package. Note: To use
      #             normalize_counts with its pseudocount parameter changed from the
      #             default 0.1 value to some positive numeric `new_var`, use:
      #             normalizer=lambda counts: doubletdetection.normalize_counts(counts,
      #             pseudocount=new_var)
      # 
      #     Attributes:
      #         all_p_values_ (ndarray): Hypergeometric test p-value per cell for cluster
      #             enrichment of synthetic doublets. Shape (n_iters, num_cells).
      #         all_scores_ (ndarray): The fraction of a cell's cluster that is
      #             synthetic doublets. Shape (n_iters, num_cells).
      #         communities_ (ndarray): Cluster ID for corresponding cell. Shape
      #             (n_iters, num_cells).
      #         labels_ (ndarray, ndims=1): 0 for singlet, 1 for detected doublet.
      #         parents_ (list of sequences of int): Parent cells' indexes for each
      #             synthetic doublet. A list wrapping the results from each run.
      #         suggested_score_cutoff_ (numeric): Cutoff used to classify cells when
      #             n_iters == 1 (scores >= cutoff). Not produced when n_iters > 1.
      #         synth_communities_ (sequence of ints): Cluster ID for corresponding
      #             synthetic doublet. Shape (n_iters, num_cells * boost_rate).
      #         top_var_genes_ (ndarray): Indices of the n_top_var_genes used. Not
      #             generated if n_top_var_genes <= 0.
      #         voting_average_ (ndarray): Fraction of iterations each cell is called a
      #             doublet.
      boost_rate <<- boost_rate
      replace <<- replace
      n_iters <<- n_iters
      normalizer <<- normalizer
      
      if(n_components == 30 & n_top_var_genes > 0){
        # If user did not change n_components, silently cap it by n_top_var_genes if needed
        n_components <<- min(n_components, n_top_var_genes)
      } else {
        n_components <<- n_components
      }
      
      # Floor negative n_top_var_genes by 0
      n_top_var_genes <<- as.integer(max(0, n_top_var_genes))
      
      #check new_lib_as function
      if(is.function(new_lib_as)){
        new_lib_as <<- new_lib_as
        print(paste("function", deparse(substitute(new_lib_as)), "accepted as new_lib_as"))
      } else if(is.null(new_lib_as) || new_lib_as == TRUE){
        new_lib_as <<- sum
        print(paste("function", 'sum', "accepted as new_lib_as"))
      } else if(new_lib_as == FALSE){
        new_lib_as <<- max
        print(paste("function", "max", "accepted as new_lib_as"))
      } else {
        stop("no valid new_lib_as input \n please enter NULL, TRUE, FALSE, or a valid function")
      }

      #check prune argument defined (to pass to phenograph)
      if(!("prune" %in% names(phenograph_parameters))){
        phenograph_parameters$prune <<- TRUE
      }
      phenograph_parameters <<- phenograph_parameters
      if(n_iters == 1 & phenograph_parameters$prune == TRUE){
        warning("Using phenograph parameter prune=FALSE is strongly recommended when \n running only one iteration. Otherwise, expect many NaN labels.")
      }
      
      if(replace == FALSE & boost_rate > 0.5){
        warning("boost_rate is trimmed to 0.5 when replace=FALSE. \n Set replace=TRUE to use greater boost rates.")
        boost_rate <<- 0.5
      }
      
      if(!((n_top_var_genes == 0) || (n_components <= n_top_var_genes))){
        stop(paste0("n_components=", n_components, " cannot be larger than n_top_var_genes=", n_top_var_genes))
      }
      print("Reference Class BoostClassifier has been initialized")
    },
    fit = function(raw_counts){
      "Fits the classifier on raw_counts."
      ##' @param raw_counts (array-like): Count matrix, oriented cells by genes.
      
      #         Args:
      #             raw_counts (array-like): Count matrix, oriented cells by genes.
      # 
      #         Sets:
      #             all_scores_, all_p_values_, communities_, top_var_genes, parents,
      #             synth_communities
      # 
      #         Returns:
      #             The fitted classifier.
      if(!is.matrix(raw_counts) & length(dim(raw_counts)) == 2){ # Only catches sparse error. Non-finite & n_dims still raised.
        if(is(raw_counts, "sparseMatrix")){
          warning("Sparse raw_counts is automatically densified")
          raw_counts <- as.matrix(raw_counts)
        } else if(is.data.frame(raw_counts)){
          warning("raw_counts data.frame automatically converted to type matrix")
          raw_counts <- as.matrix(raw_counts)
        } else {
          warning("raw_counts requires matrix input \n Attempting to convert to type matrix.")
          raw_counts <- as.matrix(raw_counts)
        }
      } else {
        if(all(is.infinite(raw_counts))){
          warning("raw_counts contains Non-finite values")
        }
      }
      
      if(n_top_var_genes > 0){
        gene_variances <- apply(raw_counts, 1, var)
        top_var_indexes <- sort(gene_variances)
        if(n_top_var_genes < nrow(raw_counts)){
          top_var_genes_ <- top_var_indexes[1:n_top_var_genes]
          #filter to top genes
          raw_counts <- raw_counts[, top_var_genes_]
        } else {
          warning("n_top_var_genes exceeds total genes \n processing full dataset")
          top_var_genes_ <- top_var_indexes[1:min(n_top_var_genes, nrow(raw_counts))]
        }
      }
      #initialise self object
      #self <- list()
      raw_counts <- raw_counts
      num_genes <- nrow(raw_counts)
      num_cells <- ncol(raw_counts)
      
      all_scores_ <- matrix(0, n_iters, num_cells)
      all_p_values_ <- matrix(0, n_iters, num_cells)
      all_communities <- matrix(0, n_iters, num_cells)
      
      all_parents <- list()
      all_synth_communities <- matrix(0, n_iters,as.integer(boost_rate * num_cells))
      
      for(i in 1:n_iters){
        print(paste0("Iteration ", i, "/", n_iters))
        all_scores_[i]  <- one_fit()
        all_p_values_[i] <- one_fit()
        all_communities[i] <- communities_
        all_parents <- c(all_parents, parents_)
        all_synth_communities[i] <- synth_communities_
      }
      # Release unneeded large data vars
      rm(raw_counts, norm_counts, rawsynthetics, synthetics)
      #attr(self, "raw_counts") <- NULL
      #attr(self, "norm_counts") <- NULL
      #attr(self, "rawsynthetics") <- NULL
      #attr(self, "synthetics") <- NULL
      
      communities_ <- all_communities
      parents_ <- all_parents
      synth_communities_ <- all_synth_communities
      
      self <- list(all_scores_, all_p_values_, communities_, top_var_genes, parents, synth_communities_)
      #return(self)
    },
    predict = function(p_thresh = 0.99, voter_thresh = 0.9){ 
      "Produce doublet calls from fitted classifier."
      #         Args:
      #             p_thresh (numeric, optional): hypergeometric test p-value threshold
      #                 that determines per iteration doublet calls
      #             voter_thresh (numeric, optional): fraction of iterations a cell must
      #                 be called a doublet
      # 
      #         Sets:
      #             labels_ and voting_average_ if n_iters > 1.
      #             labels_ and suggested_score_cutoff_ if n_iters == 1.
      # 
      #         Returns:
      #             labels_ (ndarray, ndims=1):  0 for singlet, 1 for detected doublet
      if(n_iters > 1){
        voting_average_ <- apply(all_p_values_, 1, function(x) mean(x, na.rm = TRUE) > p_thresh)
        labels_ <- ifelse(voting_average_ >= voter_thresh, voting_average_ >= voter_thresh, NA)
        voting_average_ <- ifelse(voting_average_, voting_average_, NA)
      } else{
        # Find a cutoff score
        potential_cutoffs <- unique(all_scores_[is.na(all_scores_) == FALSE])
        if(length(potential_cutoffs) > 1){
          dropoff <- potential_cutoffs[2:length(potential_cutoffs)] - potential_cutoffs[1:(1-length(potential_cutoffs))]
          max_dropoff <- which(dropoff == max(dropoff)) #+ 1 (not needed for 1-indexed language)
        } else {
          # Most likely pathological dataset, only one (or no) clusters
          max_dropoff <- 1
          suggested_score_cutoff_ <- potential_cutoffs[max_dropoff]
          labels_ <- all_scores_[1,] >= suggested_score_cutoff_ #Allow NA values
        }
      }
      return(labels_)
    },
    one_fit = function(){
      print("\nCreating downsampled doublets...")
      createDoublets()
      
      # Normalize combined augmented set
      print("Normalizing...")
      
      aug_counts <- normalizer(rbind(raw_counts, rawsynthetics_temp)) #remove _temp for internal variables
      norm_counts <- aug_counts[1:num_cells]
      synthetics <- aug_counts[num_cells:length(aug_counts)]
      
      print("Running PCA...")
      # Get phenograph results
      pca <- prcomp(aug_counts, rank = n_components, center = TRUE, scale. = TRUE)$x
      reduced_counts <- pca #apply(pca, 2, normalizer) #already normalized
      
      print("Clustering augmented data set with Phenograph...\n")
      fullcommunities <- Rphenograph(reduced_counts, k = n_components, prune = phenograph_parameters$prune)
      min_ID <- min(sizes(fullcommunities[[2]]))
      communities_ <- membership(fullcommunities[[2]])[1:num_cells]
      synth_communities_ <- membership(fullcommunities[[2]])[num_cells:length(membership(fullcommunities[[2]]))]
      community_sizes <- sizes(fullcommunities[[2]])
      
      for(ii in 1:length(community_sizes)){
        print(paste0("Found community ", names(community_sizes)[ii], " with size: ", community_sizes[ii]))
      }
      
      # Count number of fake doublets in each community and assign score
      # Number of synth/orig cells in each cluster.
      synth_cells_per_comm <- table(synth_communities_)
      orig_cells_per_comm <- table(communities_)
      community_IDs <- names(orig_cells_per_comm)
      community_scores  <- as.numeric(synth_cells_per_comm) / (synth_cells_per_comm + orig_cells_per_comm)
      scores <- sapply(1:length(communities_), function(i) community_scores[i])
      community_p_values <- sapply(1:length(community_p_values), function(i){
        phyper(synth_cells_per_comm[i], nrow(aug_counts), nrow(synthetics_temp), synth_cells_per_comm[i] + orig_cells_per_comm[i])
      })
      p_values <- sapply(1:length(communities_), function(i) community_p_values[i])
      
      if(min_ID < 0){
        scores[communities_ == -1] <- NA
        p_values[communities_ == -1] <- NA
      }
      return(list(scores, p_values))
    },
  downsampleCellPair = function(cell1, cell2){ #Downsample the sum of two cells' gene expression profiles.
    #         Args:
    #             cell1 (ndarray, ndim=1): Gene count vector.
    #             cell2 (ndarray, ndim=1): Gene count vector.
    # 
    #         Returns:
    #             ndarray, ndim=1: Downsampled gene count vector.
    
    #create doublet cells
    new_cell <- cell1 + cell2
    
    lib1 <- sum(cell1)
    lib2 <- sum(cell2)
    
    new_lib_size <- as.integer(new_lib_as(c(lib1, lib2)))
    
    mol_ind <- sample(as.integer(lib1 + lib2), size = new_lib_size)
    #mol_ind <- mol_ind #+1 #(vectorised) #not needed for 1-index language
    
    bins <- c(0, cumsum(new_cell))
    new_cell <- hist(mol_ind, bins)$counts #extract counts from histogram 
    
    return(new_cell)
  },
  createDoublets = function(){ #Create synthetic doublets.
    #         Sets .parents_
    
    # Number of synthetic doublets to add
    num_synths <- as.integer(boost_rate * num_cells)
    synthetic <- matrix(0, num_synths, num_genes)
    
    parents <- list()
    
    choices <- matrix(sample(num_cells, size=num_synths*2, replace=replace), num_synths, 2)
    synthetic <- rep(NA, nrow(choices))
    for(i in 1:nrow(choices)){
      parent_pair <- choices[i, ]
      row1 <- parent_pair[1]
      row2 <- parent_pair[2]
      if(!(is.null(new_lib_as))){
        print(paste("running dowsamplePair for synthetic cell", i))
        new_row <- downsampleCellPair(raw_counts[row1], raw_counts[row2])
      } else {
        new_row <- raw_counts[row1] + raw_counts[row2]
      }
      synthetic[i] <- new_row
      parents[[i]] <- c(row1, row2)
    }
    rawsynthetics <- synthetic
    parents_ <- parents
  }
  )
)