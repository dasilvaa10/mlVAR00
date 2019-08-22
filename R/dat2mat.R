#' @Title Function restructure a list of model coefficients into an appropriately structure matrix.
#' @name dat2mat
#' @author Alex daSilva
#' @return A pairwise matrix of coefficients
#' 
#' @param x a list of lists of model coefficients
#' @param ind how many variables in the model, used to set dimensions for output matrix
#' @param varnames names of variables, used as column and row names for output matrix
#' @param temporal are the coefficients coming in from a temporal network model?


dat2mat <- function(x, ind, temporal = FALSE, varNames) {
  
  #if not temporal (i.e., contemporaneous/between) need to tend to diagonal
  
  if (temporal == FALSE){
    
    mat <- matrix(nrow = ind, ncol = ind)
    
    diag(mat) <- 0
    
    mat <- c(mat)
    
    mat_nas <- which(is.na(mat))
    
    mat[mat_nas] <- unlist(x)
    
    mat <- matrix(mat, nrow = ind, ncol = ind, byrow = TRUE)
    
  } else {
    
    mat <- matrix(unlist(x), nrow = ind, ncol = ind, byrow = TRUE )
    
  }
  
  dimnames(mat) <- list(varNames, varNames)
  
  return(mat)
}



