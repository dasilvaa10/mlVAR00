dat2mat <- function(x, ind, temporal = FALSE, varNames) {
  
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



