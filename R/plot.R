##' @name convergence
##' @rdname plot_doublets
##'
##' @title Convergence plot of doublets
##' 
##' @description Produce a convergence plot showing number of cells called doublet per iteration
##' 
##' @param clf (BoostClassifier object): Fitted classifier
##' @param show (logical, optional): If TRUE, runs plt.show()
##' @param save (charater, optional): filename for saved figure, figure not saved by default
##' @param p_thresh (numeric, optional): hypergeometric test p-value threshold that determines per iteration doublet calls
##' @param voter_thresh (numeric, optional): fraction of iterations a cell must be called a doublet
##' 
##' @keywords scRNA quality filter matrix normalise doublets single-cell
##' @export
convergence <- function(clf, show=FALSE, save=NULL, p_thresh=0.99, voter_thresh=0.9){
    # Returns:
    #     matplotlib figure
    
    doubs_per_run <- list()
    for(i in 1:clf$n_iters){
      cum_p_values <- clf$all_p_values_[1:i]
      cum_vote_average <- apply(sum_p_values_, 2, function(x){
        x <- ifelse(is.infinite(x), NA, x)
        mean(as.numeric(x > p_thresh), na.rm = TRUE)
      })
      cum_doublets <- ifelse(cum_voting_average_ >= voter_thresh, 1, 0)
      cum_voting_average_ <- ifelse(as.numeric(cum_voting_average_), cum_voting_average_, NA)
      doubs_per_run[[i]] <- sum(cum_doublets, na.rm = TRUE)
    }

    if(show){
      plot(1:length(doubs_per_run), doubs_per_run,
           type = "l", col = "royalblue2",
           main = "Predicted Doublets\n per Iteration",
           xlab = "Number of Predicted Doublets",
           ylab = "Number of Iterations")
    }
    
    if(is.character(save)){
      if(strsplit(save, split = "[.]")[[1]][2] == "png"){
        save <- save
        print(paste("Saving tSNE plot as file:", save))
        png(file = save, width = 800, height = 800)
          plot(1:length(doubs_per_run), doubs_per_run,
              type = "l", col = "royalblue2",
              main = "Predicted Doublets\n per Iteration",
              xlab = "Number of Predicted Doublets",
              ylab = "Number of Iterations")
        dev.off()
      } else {
        if(strsplit(save, split = "[.]")[[1]][2] == "pdf"){
          save <- save
          print(paste("Saving tSNE plot as file:", save))
        } else {
          save <- strsplit(save, split = "[.]")[[1]]
          save <- save[1:(length(save)-1)] # remove extension
          save <- paste0(save, ".pdf")
          warning("file extension in save changed to pdf")
          print(paste("Saving tSNE plot as file:", save))
        }
        pdf(file = save, width = 8, height = 8)
          plot(1:length(doubs_per_run), doubs_per_run,
              type = "l", col = "royalblue2",
              main = "Predicted Doublets\n per Iteration",
              xlab = "Number of Predicted Doublets",
              ylab = "Number of Iterations")
        dev.off()
      }
    } else {
      warning("Convergence plot not saved, give a valid filename for save")
    }
}
##' @rdname plot_doublets
##' 
##' @title  tSNE Plot
##' 
##' @description Produce a tsne plot of the data with doublets in black
##' 
##' @param raw_counts (matrix): genes by cells count matrix
##' @param labels (vector, integer): predicted doublets from predict method
##' @param n_components (integer, optional): number of PCs to use prior to TSNE
##' @param n_jobs (integer, optional): number of cores to use for TSNE, -1 for all
##' @param show (logical, optional): If TRUE, runs plt.show()
##' @param save (charater, optional): filename for saved figure, figure not saved by default
##' @import gplots Rtsne Rcpp
##' 
##' @export
tsne_plot <- function(raw_counts, labels, n_components=30L, n_jobs=-1, show=False, save=None){
      #     """Produce a tsne plot of the data with doublets in black
      #
      #     Returns:
      #         matplotlib figure
      #         ndarray: tsne reduction
      #     """
      norm_counts <- normalize_counts(raw_counts)
      
      pca <- prcomp(aug_counts, rank = n_components, center = TRUE, scale. = TRUE)$rotation
      reduced_counts <- pca #apply(pca, 2, normalizer) #already normalized
      
      communities <- Rphenograph(reduced_counts, k = n_components, prune = phenograph_parameters$prune)[[2]]
      
      #reduced_counts <- fit_transform(reduced_counts)                     
      
      tsne_counts <- Rtsne(reduced_counts)
      
      plot(tsne[,1], tsne[,2], col = ifelse(labels, rainbow(n_components)[communities], "black"),
           main = "Cells with Detected\n Doublets in Black",
           sub = paste(sum(labels), "doublets out of", ncol(raw_counts),  "cells.\n",  round(100 * nsum(labels) / nrow(raw_counts.shape), 2),  "across-type doublet rate"),
           xaxt = "n", yaxt = "n")
  
    
      if(show){
        plot(tsne[,1], tsne[,2], col = ifelse(labels, rainbow(n_components)[communities], "black"),
             main = "Cells with Detected\n Doublets in Black",
             sub = paste(sum(labels), "doublets out of", ncol(raw_counts),  "cells.\n",  round(100 * nsum(labels) / nrow(raw_counts.shape), 2),  "across-type doublet rate"),
             xaxt = "n", yaxt = "n")
      }
      if(is.character(save)){
        if(strsplit(save, split = "[.]")[[1]][2] == "png"){
          save <- save
          print(paste("Saving tSNE plot as file:", save))
          png(file = save, width = 800, height = 800)
            plot(tsne[,1], tsne[,2], col = ifelse(labels, rainbow(n_components)[communities], "black"),
                 main = "Cells with Detected\n Doublets in Black",
                 sub = paste(sum(labels), "doublets out of", ncol(raw_counts),  "cells.\n",  round(100 * nsum(labels) / nrow(raw_counts.shape), 2),  "across-type doublet rate"),
                 xaxt = "n", yaxt = "n") 
          dev.off()
        } else {
          if(strsplit(save, split = "[.]")[[1]][2] == "pdf"){
            save <- save
            print(paste("Saving tSNE plot as file:", save))
          } else {
            save <- strsplit(save, split = "[.]")[[1]]
            save <- save[1:(length(save)-1)] # remove extension
            save <- paste0(save, ".pdf")
            warning("file extension in save changed to pdf")
            print(paste("Saving tSNE plot as file:", save))
          }
          pdf(file = save, width = 8, height = 8)
          plot(tsne[,1], tsne[,2], col = ifelse(labels, rainbow(n_components)[communities], "black"),
               main = "Cells with Detected\n Doublets in Black",
               sub = paste(sum(labels), "doublets out of", ncol(raw_counts),  "cells.\n",  round(100 * nsum(labels) / nrow(raw_counts.shape), 2),  "across-type doublet rate"),
               xaxt = "n", yaxt = "n") 
          dev.off()
        }
      } else {
        warning("tSNE plot not saved, give a valid filename for save")
      }
      return(tsne_counts)
}