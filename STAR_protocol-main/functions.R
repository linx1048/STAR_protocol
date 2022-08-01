#' Depth log-normalization
#'
#' 1) Normalizes reads by sequencing depth with a pseudocount and scaling factor
#' 2) Log2-normalization  
#' 
#' @param df Data frame of read counts.
#' @param scale Scaling factor (default 1 million).
#' @param pseudo Pseudocount (default 1).
#' @return Log- and depth-normalized read counts. (NA removed during normalization)
#' @export
normalize_reads <- function(df, scale = 1e6, pseudo = 1) {
  log2((df / sum(df, na.rm = TRUE)) * scale + pseudo)
}

#' LFC
#'
#' Calculate average T0, then calculate LFC relative to avg. T0
#' 
#' @param df Data frame of normalized read counts.
#' @param t0 numerical column denoting T0 samples.
#' @return guide-level LFC values 
#' @export
calc_lfc <- function(df, t0 = 1:3) {
  df1 = df[,-c(1:2)]
  avg_t0 = rowMeans(df1[,t0])
  df2 = df1[,-t0]
  df3 = df2 - avg_t0
  df4 = cbind(df[,c(1:2)], df3)
  return(df4)
  
}

#' Computes essential gene recovery AUC.
#' 
#' Computes area under the curve for ROC curves that measure how well each technical replicate
#' recovers signal for essential-targeting guides.
#' 
#' @param df LFC dataframe.
#' @return Returns a dataframe with three columns for replicate name, essential AUC relative 
#'   to all other genes
essential_qc <- function(df, essentials_set, screens) {
  
  # Checks that the given gene_col is in the data
  if (!("gene" %in% colnames(df))) {
    stop(paste("gene name column 'gene' not in df"))
  }
  
  # Loads essential gene standard from internal data
  essentials <- essentials_set
  
  # Gets indices of essential-targeting guides
  essential_ind <- df$gene %in% essentials
  
  # Throws warning if too few genes in standards
  if (sum(essential_ind) < 10) {
    warning(paste("too few essential-targeting guides in df, skipping all AUC computation"))
    return(NULL)
  }
  
  # Gets PR curves for all essential genes and all technical replicates
  results <- data.frame(screen = NA, 
                        essential_AUC_all = NA)
  counter <- 1
  for (rep in screens){
      rep_name = colnames(df)[rep]
      temp <- df[,c(1,rep)]
      ind <- temp$gene %in% essentials
      if (sum(ind) < 10) {
        warning(paste("too few essential-targeting guides for replicate", rep, ", skipping AUC computation"))
      }
      precrec_obj <- evalmod(scores = -temp[[2]], labels = as.numeric(ind))
      autoplot(precrec_obj, curvetype = "ROC", show_legend = F)
      auc1 <- attributes(precrec_obj)$auc$aucs[1]
      results[counter,] <- c(rep_name, auc1)
      counter <- counter + 1


      
    }

  return(results)
}
