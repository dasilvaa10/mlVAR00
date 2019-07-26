#' @Title Function to force a matrix to be positive.
#' @name forcePositive
#' @author Sacha Epskamp,  https://github.com/SachaEpskamp/mlVAR
#' @return A positive matrix
#' 
#' @param x data matrix


forcePositive <- function(x){
  x <- (x + t(x))/2
  
  if (any(eigen(x)$values < 0)){
    
    return(x - (diag(nrow(x)) * min(eigen(x)$values)-0.001)) 
    
  } else {
    
    return(x)
    
  }
  
}
