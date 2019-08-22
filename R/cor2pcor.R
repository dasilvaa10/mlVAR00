#' @Title Function to convert a correlation matrix into a partial correlation matrix.
#' @name cor2pcor
#' @author Alex daSilva
#' @return A matrix containing partial correlations
#' 
#' @param x either raw data or a correlation matrix
#' @param type indicates whether the data are raw or a correlation matrix
#' @importFrom MASS ginv
 

PCOR <- function(x, type = c("raw", "cor")) {
  
  #check type, if raw, calculate correlation matrix
  
  type <- match.arg(type)
  
  if (type == "raw") {
    
    x <- scale(x)
    
    R <- (t(x) %*% x) / (nrow(x) - 1)
    
  } else  {
    
    R <- x
    
  }
  
  #calulcate anti-image correlation matrix
  
  ind <- unique(dim(R))
  
  R_inv <- ginv(R)
  
  ZM <- matrix(rep(0, len = (ind*ind)), nrow = ind)
  
  diag(ZM) <- diag(R_inv)
  
  D <- ginv(ZM)
  
  AICOV <- D %*% R_inv %*% D
  
  diag(ZM) <- diag(AICOV)
  
  D  <- ginv(sqrt(ZM))
  
  AICOR <- D %*% AICOV %*% D
  
  pcor <- AICOR
  
  #multiply a negative through to get partial correlations
  
  pcor[upper.tri(pcor)] <- -pcor[upper.tri(pcor)]
  
  pcor[lower.tri(pcor)] <- -pcor[lower.tri(pcor)]
  
  dimnames(pcor) <- list(colnames(R), colnames(R))
  
  return(pcor)
  
}  
