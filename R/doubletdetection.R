##' @name normalize_counts
##' @rdname normalize_counts
##'
##' @title Normalize count array.
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
##' @description Functions required to compute identification of doublets in single-cell RNA-Seq experiments.
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
  raw_counts <<- readMM(paste(group, "matrix.mtx", sep = "/"))
  rownames(raw_counts) <- gene_names
  colnames(raw_counts) <- barcodes
  raw_counts <<- as.matrix(raw_counts)
  return(raw_counts)
}

##' @rdname load_10x_h5
##' 
##' @title  Load count matrix in 10x H5 format
##' 
##' @description Functions required to compute identification of doublets in single-cell RNA-Seq experiments.
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
         "'. Please double-check if the directory exists.")
  if (!dir.exists(normalizePath(paste(dirname(file))))) 
    stop("Could not find the pipestance output directory: '", 
         normalizePath(paste(dirname(file))), "'. Please double-check if the directory exists.")
  h5_path <- file.path(normalizePath(paste(dirname(file))), ifelse(barcode_filtered, 
                                                                   "filtered_gene_bc_matrices_h5.h5", "raw_gene_bc_matrices_h5.h5"))
  if (!file.exists(h5_path)) {
    stop(paste("Could not find matrix H5:", h5_path))
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
##' @import Matrix stats Rphenograph igraph
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
##' @field num_genes,num_cells (numeric): number of genes and cells, rows and columns of raw_counts matrix respectively.
##' @field rawsynthetics (matrix): Count matrix, oriented genes by cells. Synthetic doublets generated.
##' @field parents_ (list): Pairs of column indices for cells sampled to generate synthetic doublets.
##' 
##' 
##' @param raw_counts (matrix): Count matrix, oriented genes by cells. Can be sparse matrix or data.frame input.
##' @param p_thresh (numeric, optional): hypergeometric test p-value threshold that determines per iteration doublet calls
##' @param voter_thresh (numeric, optional): fraction of iterations a cell must be called a doublet
##' @param cell1,cell2 (vector, numeric): Gene count vectors.

##' @usage
##' 
##' library("DoubletDetection")
##' clf <- BoostClassifier$new()
##' # raw_counts is a cells by genes count matrix
##' labels <- clf$fit(raw_counts)$predict()
##' #returns a vector of 1 for doublet and 0 for singlet
##' 
##' @export BoostClassifier
##' @exportClass BoostClassifier

BoostClassifier <- setRefClass(
  "BoostClassifier",
  fields = list(
    boost_rate = "numeric", #ANY for no type specified (all passed on)
    n_components = "numeric",
    n_top_var_genes = "numeric",
    new_lib_as = "ANY",
    replace = "logical",
    phenograph_parameters = "list",
    n_iters = "numeric",
    normalizer = "function",
    num_genes = "numeric",
    num_cells = "numeric",
    rawsynthetics = "matrix",
    all_scores_ = "ANY", 
    all_p_values_ = "ANY", #DEPRECATED
    all_log_p_values_ = "ANY",
    communities_ = "ANY",
    top_var_genes_ = "ANY",
    parents_ = "list",
    synth_communities_ = "ANY",
    raw_counts = "ANY"#,
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
      #             enrichment of synthetic doublets. Shape (n_iters, num_cells). Due to rounding point errors,
      #             use of all_log_p_values recommended. Will be removed in v3.0.
      #         all_log_p_values_ (ndarray): Hypergeometric test natural log p-value per
      #            cell for cluster enrichment of synthetic doublets. Shape (n_iters,num_cells).
      #         all_scores_ (ndarray): The fraction of a cell's cluster that is
      #             synthetic doublets. Shape (n_iters, num_cells).
      #         communities (ndarray): Cluster ID for corresponding cell. Shape
      #             (n_iters, num_cells).
      #         labels_ (ndarray, ndims=1): 0 for singlet, 1 for detected doublet.
      #         parents_ (list of sequences of int): Parent cells' indexes for each
      #             synthetic doublet. A list wrapping the results from each run.
      #         suggested_score_cutoff_ (numeric): Cutoff used to classify cells when
      #             n_iters == 1 (scores >= cutoff). Not produced when n_iters > 1.
      #         synth_communities (sequence of ints): Cluster ID for corresponding
      #             synthetic doublet. Shape (n_iters, num_cells * boost_rate).
      #         top_var_genes_ (ndarray): Indices of the n_top_var_genes used. Not
      #             generated if n_top_var_genes <= 0.
      #         voting_average_ (ndarray): Fraction of iterations each cell is called a
      #             doublet.
      if(!is.integer(n_components)){
        n_components <<- as.integer(n_components)
        if(n_components == floor(n_components)){
          warning("numeric input for n_components taken as an integer")
        } else{
          n_components <<- floor(n_components)
          warning("numeric input for n_components rounded down an integer")
        }
      }
      if(!is.integer(n_top_var_genes)){
        n_top_var_genes <<- as.integer(n_top_var_genes)
        if(n_top_var_genes == floor(n_top_var_genes)){
          warning("numeric input for n_top_var_genes taken as an integer")
        } else{
          n_top_var_genes <<- floor(n_top_var_genes)
          warning("numeric input for n_top_var_genes rounded down an integer")
        }
      }
      if(!is.integer(n_iters)){
        n_iters <<- as.integer(n_iters)
        if(n_iters == floor(n_iters)){
          warning("numeric input for n_iters taken as an integer")
        } else{
          n_iters <<- floor(n_iters)
          warning("numeric input for n_iters rounded down an integer")
        }
      }
      
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
      
      # Floor negative n_top_var_genes by 1 (1-indexed language)
      n_top_var_genes <<- as.integer(max(1, n_top_var_genes))
      
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
        stop("no valid new_lib_as input
              please enter NULL, TRUE, FALSE, or a valid function")
      }

      #check prune argument defined (to pass to phenograph)
      if(!("prune" %in% names(phenograph_parameters))){
        phenograph_parameters$prune <<- TRUE
      }
      phenograph_parameters <<- phenograph_parameters
      if(n_iters == 1 & phenograph_parameters$prune == TRUE){
        warning("Using phenograph parameter prune=FALSE is strongly recommended when
                running only one iteration. Otherwise, expect many NaN labels.")
      }
      
      if(replace == FALSE & boost_rate > 0.5){
        warning("boost_rate is trimmed to 0.5 when replace=FALSE.
                Set replace=TRUE to use greater boost rates.")
        boost_rate <<- 0.5
      }
      
      if(!((n_top_var_genes == 0) || (n_components <= n_top_var_genes))){
        stop(paste0("n_components=", n_components, " cannot be larger than n_top_var_genes=", n_top_var_genes))
      }
      print("Reference Class BoostClassifier has been initialized")
    },
    fit = function(raw_counts){
      "Fits the classifier on raw_counts."
      #         Args:
      #             raw_counts (array-like): Count matrix, oriented cells by genes.
      # 
      #         Sets:
      #             all_scores_, all_p_values_, all_log_p_values, communities, top_var_genes, parents,
      #             synth_communities
      # 
      #         Returns:
      #             The fitted classifier.
      if(!is.matrix(raw_counts) & length(dim(raw_counts)) == 2){ # Only catches sparse error. Non-finite & n_dims still raised.
        if(is(raw_counts, "sparseMatrix")){
          warning("Sparse raw_counts is automatically densified")
          raw_counts <<- as.matrix(raw_counts)
        } else if(is.data.frame(raw_counts)){
          warning("raw_counts data.frame automatically converted to type matrix")
          raw_counts <<- as.matrix(raw_counts)
        } else {
          warning("raw_counts requires matrix input")
          print("Attempting to convert to type matrix.")
          raw_counts <<- as.matrix(raw_counts)
        }
      } else {
        if(all(is.infinite(raw_counts))){
          warning("raw_counts contains Non-finite values")
        }
      }
      
      if(n_top_var_genes > 0){
        gene_variances <- apply(raw_counts, 1, var)
        top_var_indexes <- order(gene_variances)
        if(n_top_var_genes < nrow(raw_counts)){
          top_var_genes_ <<- top_var_indexes[1:n_top_var_genes]
          #filter to top genes
          raw_counts <<- raw_counts[top_var_genes_,]
        } else {
          warning("n_top_var_genes exceeds total genes")
          print("processing full dataset")
          top_var_genes_ <<- top_var_indexes[1:min(n_top_var_genes, nrow(raw_counts))]
        }
      }
      #initialise self object
      #self <- list()
      raw_counts <<- raw_counts
      num_genes <<- nrow(raw_counts)
      num_cells <<- ncol(raw_counts)
      
      all_scores_ <<- matrix(0, n_iters, num_cells)
      all_p_values_ <<- matrix(0, n_iters, num_cells)
      all_log_p_values_ <<- matrix(0, n_iters, num_cells)
      all_communities <- matrix(0, n_iters, num_cells)
      
      all_parents <- list()
      all_synth_communities <- matrix(0, n_iters, as.integer(boost_rate * num_cells))
      
      for(i in 1:n_iters){
        print(paste0("Iteration ", i, "/", n_iters))
        fits <- one_fit()
        all_scores_[i,]  <<- fits$scores
        all_p_values_[i,] <<- fits$p_values #DEPRECATED
        all_log_p_values_[i,] <<- fits$log_p_values
        all_communities[i,] <- fits$communities
        all_parents[[i]] <- fits$parents_
        all_synth_communities[i,] <- fits$synth_communities
      }
      # Release unneeded large data vars
      #rm(raw_counts, norm_counts, rawsynthetics, synthetics)
      #attr(self, "raw_counts") <- NULL
      #attr(self, "norm_counts") <- NULL
      #attr(self, "rawsynthetics") <- NULL
      #attr(self, "synthetics") <- NULL
      
      communities_ <<- all_communities
      parents_ <<- all_parents
      synth_communities_ <<- all_synth_communities
      
      self <- list(all_scores_, all_p_values_, all_log_p_values_, communities_, top_var_genes_, parents_, synth_communities_)
      names(self) <- c("all_scores_", "all_p_values_", "all_log_p_values_", "communities_", "top_var_genes_", "parents_", "synth_communities_")
      return(.self)
    },
    predict = function(p_thresh = 0.01, voter_thresh = 0.9){ 
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
      log_p_thresh <- log(p_thresh)
      if(n_iters > 1){
        voting_average_ <- apply(all_log_p_values_, 2, function(x){
          x <- ifelse(is.infinite(x), NA, x)
          mean(as.numeric(x <= log_p_thresh), na.rm = TRUE)
        })
        labels_ <- ifelse(voting_average_ >= voter_thresh, 1, 0)
        voting_average_ <- ifelse(is.numeric(voting_average_), voting_average_, NA)
      } else{
        # Find a cutoff score
        potential_cutoffs <- unique(sort(all_scores_[is.na(all_scores_) == FALSE], decreasing = TRUE))
        if(length(potential_cutoffs) > 1){
          dropoff <- potential_cutoffs[2:length(potential_cutoffs)] - potential_cutoffs[1:(length(potential_cutoffs)-1)]
          max_dropoff <- which(dropoff == max(dropoff)) #+ 1 (not needed for 1-indexed language)
        } else {
          # Most likely pathological dataset, only one (or no) clusters
          max_dropoff <- 1
        }
        suggested_score_cutoff_ <- potential_cutoffs[max_dropoff]
        labels_ <- ifelse(all_scores_ >= suggested_score_cutoff_, 1, 0) #Allow NA values
      }
      return(labels_)
    },
    one_fit = function(){
      print("Creating downsampled doublets...")
      createDoublets()
      
      # Normalize combined augmented set
      print("Normalizing...")
      
      #print(dim(raw_counts))
      #print(dim(rawsynthetics))
      
      aug_counts <- cbind(raw_counts, rawsynthetics) #combine raw_counts and synthetic doublets
      gene_counts <- apply(aug_counts, 1, sum)
      if(any(gene_counts == 0)){
        aug_counts <- aug_counts[apply(aug_counts, 1, sum) != 0,] #remove genes with zero total counts
        warning(paste(sum(gene_counts == 0), "genes with zero total removed:", sum(gene_counts != 0), "genes remaining"))
      }
      aug_counts <- normalizer(aug_counts) #normalise combined ataset
      
      norm_counts <- aug_counts[,1:num_cells]
      synthetics <- aug_counts[, (num_cells + 1):ncol(aug_counts)]
      
      print("Running PCA...")
      # Get phenograph results
      pca <- prcomp(t(aug_counts), rank = n_components, center = TRUE, scale. = TRUE)$x
      #pca <- normalizer(t(pca))
      reduced_counts <- pca[, 1:min(n_components, ncol(pca))]
      
      print("Clustering augmented data set with Phenograph...")
      fullcommunities <- Rphenograph(reduced_counts, k = n_components, prune = phenograph_parameters$prune)
      min_ID <- min(sizes(fullcommunities[[2]]))
      communities <- membership(fullcommunities[[2]])[1:num_cells]
      synth_communities <- membership(fullcommunities[[2]])[(num_cells + 1):length(membership(fullcommunities[[2]]))]
      community_sizes <- sizes(fullcommunities[[2]])
      
      for(ii in 1:length(community_sizes)){
        print(paste0("Found community ", names(community_sizes)[ii], " with size: ", community_sizes[ii]))
      }
      
      # Count number of fake doublets in each community and assign score
      # Number of synth/orig cells in each cluster.
      orig_cells_per_comm <- table(communities)
      community_IDs <- names(orig_cells_per_comm)
      synth_cells_per_comm <- table(synth_communities)
      synth_cells_per_comm <- as.table(ifelse(community_IDs %in% names(synth_cells_per_comm), synth_cells_per_comm[match(community_IDs, names(synth_cells_per_comm))], 0))
      names(synth_cells_per_comm)<- community_IDs 
      community_scores  <- as.numeric(synth_cells_per_comm) / (synth_cells_per_comm + orig_cells_per_comm)
      scores <- community_scores[communities]
      #community_p_values <- sapply(1:length(community_IDs), function(i){
      #  phyper(synth_cells_per_comm[i], ncol(synthetics), ncol(raw_counts), synth_cells_per_comm[i] + orig_cells_per_comm[i], lower.tail = FALSE)
      #})  #DEPRECATED
      community_log_p_values <- sapply(1:length(community_IDs), function(i){
        phyper(synth_cells_per_comm[i], ncol(synthetics), ncol(raw_counts), synth_cells_per_comm[i] + orig_cells_per_comm[i], lower.tail = FALSE, log.p = TRUE)
      })
      community_p_values <- exp(community_log_p_values)  #DEPRECATED
      p_values <- community_p_values[communities] #DEPRECATED
      log_p_values <- community_log_p_values[communities]
    
      if(min_ID < 0){
        scores[communities == -1] <- NA
        p_values[communities == -1] <- NA #DEPRECATED
        log_p_values[communities == -1] <- NA
      }
      
      outs <- list(scores, p_values, log_p_values, communities, synth_communities, parents_)
      names(outs) <- c("scores", "p_values", "log_p_values", "communities", "synth_communities", "parents_")
      return(outs)
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
    new_cell <- hist(mol_ind, bins, plot = FALSE)$counts #extract counts from histogram 
    
    return(new_cell)
  },
  createDoublets = function(){ #Create synthetic doublets.
    #         Sets .parents_
    
    num_genes <<- nrow(raw_counts)
    num_cells <<- ncol(raw_counts)
    #print(dim(raw_counts))
    
    # Number of synthetic doublets to add
    num_synths <- as.integer(boost_rate * num_cells)
    synthetic <- matrix(0, num_genes, num_synths)
    
    parents <- as.list(rep(NA, num_synths))
    
    choices <- matrix(sample(num_cells, size=num_synths*2, replace=replace), 2, num_synths)
    #synthetic <- matrix(NA, num_genes, ncol(choices))
    for(i in 1:ncol(choices)){
      parent_pair <- choices[, i]
      row1 <- parent_pair[1]
      row2 <- parent_pair[2]
      if(!(is.null(new_lib_as))){
        print(paste("running downsamplePair for synthetic cell", i))
        new_row <- downsampleCellPair(raw_counts[,row1], raw_counts[,row2])
      } else {
        new_row <- raw_counts[,row1] + raw_counts[,row2]
      }
      synthetic[,i] <- new_row
      parents[[i]] <- c(row1, row2)
    }
    rawsynthetics <<- synthetic
    #print(dim(synthetic))
    parents_ <<- parents
    return(list(parents_, rawsynthetics))
  }
  )
)
