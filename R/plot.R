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
convergence <- function(clf, show=FALSE, save=NULL, p_thresh=0.01, voter_thresh=0.9){
    # Returns:
    #     matplotlib figure
    
    log_p_thres <- log(p_thresh)
    doubs_per_run <- list()
    for(i in 1:clf$n_iters){
      cum_log_p_values_ <- matrix(clf$all_log_p_values_[1:i,], i, ncol(clf$all_log_p_values_))
      cum_voting_average_ <- apply(cum_log_p_values_, 2, function(x) {
        x <- ifelse(is.infinite(x), NA, x)
        mean(as.numeric(x <= log_p_thresh), na.rm = TRUE)
      })
      cum_doublets <- ifelse(cum_voting_average_ >= voter_thresh, 1, 0)
      cum_voting_average_ <- ifelse(is.numeric(cum_voting_average_), cum_voting_average_, NA)
      doubs_per_run[[i]] <- sum(cum_doublets, na.rm = T)
    }

    if(show){
      plot(1:length(doubs_per_run), doubs_per_run,
           type = "l", col = "royalblue2",
           main = "Predicted Doublets\n per Iteration",
           xlab = "Number of Iterations", 
           ylab = "Number of Predicted Doublets")
           
    }
    
    if(is.character(save)){
      if(strsplit(save, split = "[.]")[[1]][2] == "png"){
        save <- save
        print(paste("Saving tSNE plot as file:", save))
        png(file = save, width = 800, height = 800)
          plot(1:length(doubs_per_run), doubs_per_run,
              type = "l", col = "royalblue2",
              main = "Predicted Doublets\n per Iteration",
              xlab = "Number of Iterations", 
              ylab = "Number of Predicted Doublets")
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
              xlab = "Number of Iterations", 
              ylab = "Number of Predicted Doublets")
        dev.off()
      }
    } else {
      warning("Convergence plot not saved, give a valid filename for save")
    }
    return(doubs_per_run)
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
##' @importFrom graphics mtext plot title
##' @importFrom grDevices dev.off pdf png rainbow
##' 
##' @export
tsne_plot <- function(raw_counts, labels, n_components=30L, n_jobs=-1, show=FALSE, save=NULL){
      #     """Produce a tsne plot of the data with doublets in black
      #
      #     Returns:
      #         matplotlib figure
      #         ndarray: tsne reduction
      #     """
      norm_counts <- normalize_counts(raw_counts)
      
      pca <- prcomp(norm_counts, rank = n_components, center = TRUE, scale. = TRUE)$rotation
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

##' @rdname plot_doublets
##' 
##' @title  Threshold Plot
##' 
##' @description Produce a plot showing number of cells called doublet across various thresholds
##' 
##' @param clf (BoostClassifier object): Fitted classifier
##' @param log10 (logical, optional): Use natural log p values if FALSE, log10 otherwise.
##' @param n_components (integer, optional): number of PCs to use prior to TSNE
##' @param show (logical, optional): If TRUE, runs plt.show()
##' @param save (charater, optional): filename for saved figure, figure not saved by default
##' @param p_grid (numeric vector, optional): p-value thresholds to use
##' @param voter_grid (numeric vector, optional):  Voting thresholds to use. Defaults to c(0.3, 1.0, 0.01).
##' @param p_step (integer, optional): number of xlabels to skip in plot
##' @param v_step (integer, optional): number of ylabels to skip in plot
##' 
##' @import gplots heatmap.2x
##' 
##' @export

threshold <- function(clf, log10=True, show=FALSE, save=NULL, p_grid=NULL, voter_grid=NULL, v_step=20,
              p_step=20){
          # Produce a plot showing number of cells called doublet across
          # various thresholds
          # Args:
          # clf (BoostClassifier object): Fitted classifier
          # log10 (bool, optional): Use natural log p values if False, log10
          # otherwise.
          # show (bool, optional): If True, runs plt.show()
          # save (str, optional): If provided, the figure is saved to this
          # filepath.
          # p_grid (ndarray, optional): p-value thresholds to use
          # voter_grid (ndarray, optional): Voting thresholds to use. Defaults to
          # np.arange(0.3, 1.0, 0.01).
          # p_step (int, optional): number of xlabels to skip in plot
          # v_step (int, optional): number of ylabels to skip in plot
          # Returns:
          # matplotlib figure
          
  all_p_values <- clf$all_log_p_values_
  if(log10){
    log_function <- log10
    all_p_values <- all_p_values / log(10)
  } else{
    log_function <- log
  }
  if(is.null(p_grid)){
    p_grid <- unique(sort(all_p_values_))
    p_grid <- p_grid[p_grid < log_function(0.01)]
  }
  if(is.null(voter_grid)){
    voter_grid <- c(0.3, 1.0, 0.01)
  }
  voter_grid <- sort(voter_grid, decreasing = F)
  doubs_per_t <- matrix(0, length(voter_grid), length(p_grid))
  for(i in 1:length(voter_grid)){
    for(j in 1:length(p_grid)){
      voting_average <- apply(all_p_values_ <= p_grid[j], 1, mean, na.rm=T)
      labels = as.integer(voting_average >= voter_grid[i])
      doubs_per_t[i, j] = sum(labels, na.rm =T)
    }
  }
  heatmap.2x(doubs_per_t, scale = 'none', trace = 'none', 
            dendrogram = 'none', Rowv = FALSE, Colv = FALSE,
            key.title = "Predicted Doublets",
            col = colorpanel(50, "black", "red", "yellow"),
            labCol = ifelse(1:length(p_grid) %in% seq(1, length(p_grid), p_step), round(p_grid, 2), ""),
            labRow = ifelse(1:length(voter_grid) %in% seq(1, length(voter_grid), v_step), voter_grid, ""),
            cexLab = 7)
  #axis(1, at = 0.21 + seq(0, 1, p_step/(length(p_grid)-1)) * 0.72, labels = round(p_grid[seq(1, length(p_grid), p_step)], 2), line = 1, las = 2 )
  #axis(2, at = -0.05 + seq(0, 1, v_step/(length(voter_grid)-1)) * 0.85, labels = voter_grid, line = -6.5, las = 2 )
  
  if(log10){
    title(xlab = "Log10 p-value")
  } 
  else {
    title(xlab = "Log p-value")
  }
  mtext("Voting Threshold", side = 4)
  title('Threshold Diagnostics')

  if(show){
    heatmap.2x(doubs_per_t, scale = 'none', trace = 'none', 
               dendrogram = 'none', Rowv = FALSE, Colv = FALSE,
               key.title = "Predicted Doublets",
               col = colorpanel(50, "black", "red", "yellow"),
               labCol = ifelse(1:length(p_grid) %in% seq(1, length(p_grid), p_step), round(p_grid, 2), ""),
               labRow = ifelse(1:length(voter_grid) %in% seq(1, length(voter_grid), v_step), voter_grid, ""),
               cexLab = 7)
    #axis(1, at = 0.21 + seq(0, 1, p_step/(length(p_grid)-1)) * 0.72, labels = round(p_grid[seq(1, length(p_grid), p_step)], 2), line = 1, las = 2 )
    #axis(2, at = -0.05 + seq(0, 1, v_step/(length(voter_grid)-1)) * 0.85, labels = voter_grid, line = -6.5, las = 2 )
    
    if(log10){
      title(xlab = "Log10 p-value")
    } 
    else {
      title(xlab = "Log p-value")
    }
    mtext("Voting Threshold", side = 4)
    title('Threshold Diagnostics')
    
  }
  if(is.character(save)){
    if(strsplit(save, split = "[.]")[[1]][2] == "png"){
      save <- save
      print(paste("Saving threshold plot as file:", save))
      png(file = save, width = 800, height = 800)
      heatmap.2x(doubs_per_t, scale = 'none', trace = 'none', 
                 dendrogram = 'none', Rowv = FALSE, Colv = FALSE,
                 key.title = "Predicted Doublets",
                 col = colorpanel(50, "black", "red", "yellow"),
                 labCol = ifelse(1:length(p_grid) %in% seq(1, length(p_grid), p_step), round(p_grid, 2), ""),
                 labRow = ifelse(1:length(voter_grid) %in% seq(1, length(voter_grid), v_step), voter_grid, ""),
                 cexLab = 7)
      #axis(1, at = 0.21 + seq(0, 1, p_step/(length(p_grid)-1)) * 0.72, labels = round(p_grid[seq(1, length(p_grid), p_step)], 2), line = 1, las = 2 )
      #axis(2, at = -0.05 + seq(0, 1, v_step/(length(voter_grid)-1)) * 0.85, labels = voter_grid, line = -6.5, las = 2 )
      
      if(log10){
        title(xlab = "Log10 p-value")
      } 
      else {
        title(xlab = "Log p-value")
      }
      mtext("Voting Threshold", side = 4)
      title('Threshold Diagnostics')
      
      dev.off()
    } else {
      if(strsplit(save, split = "[.]")[[1]][2] == "pdf"){
        save <- save
        print(paste("Saving threshold plot as file:", save))
      } else {
        save <- strsplit(save, split = "[.]")[[1]]
        save <- save[1:(length(save)-1)] # remove extension
        save <- paste0(save, ".pdf")
        warning("file extension in save changed to pdf")
        print(paste("Saving tSNE plot as file:", save))
      }
      pdf(file = save, width = 8, height = 8)
      heatmap.2x(doubs_per_t, scale = 'none', trace = 'none', 
                 dendrogram = 'none', Rowv = FALSE, Colv = FALSE,
                 key.title = "Predicted Doublets",
                 col = colorpanel(50, "black", "red", "yellow"),
                 labCol = ifelse(1:length(p_grid) %in% seq(1, length(p_grid), p_step), round(p_grid, 2), ""),
                 labRow = ifelse(1:length(voter_grid) %in% seq(1, length(voter_grid), v_step), voter_grid, ""),
                 cexLab = 7)
      #axis(1, at = 0.21 + seq(0, 1, p_step/(length(p_grid)-1)) * 0.72, labels = round(p_grid[seq(1, length(p_grid), p_step)], 2), line = 1, las = 2 )
      #axis(2, at = -0.05 + seq(0, 1, v_step/(length(voter_grid)-1)) * 0.85, labels = voter_grid, line = -6.5, las = 2 )
      
      if(log10){
        title(xlab = "Log10 p-value")
      } 
      else {
        title(xlab = "Log p-value")
      }
      mtext("Voting Threshold", side = 4)
      title('Threshold Diagnostics')
      
      dev.off()
    }
  } else {
    warning("threshold plot not saved, give a valid filename for save")
  }
  return(doubs_per_t)
}