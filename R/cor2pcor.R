#' @Title Function to convert a correlation matrix into a partial correlation matrix.
#' @name cor2pcor
#' @author Alex daSilva, BAsed on G. Jay Kerns, Ph.D., Youngstown State University (http://tolstoy.newcastle.edu.au/R/e2/help/07/08/22816.html)
#' @return A matrix containing partial correlations
#' 
#' @param cormat a matrix of correlations
#' @importFrom MASS ginv
 


cor2pcor <- function(cormat) {
  
  X <- cormat 
  
  iX <- MASS::ginv(X) 
  
  S2 <- diag(diag((iX^-1)))
  
  # anti-image covariance matrix
  
  AIS <- S2%*%iX%*%S2 
  
  # image covariance matrix
  
  IS <- X+AIS-2*S2 
  
  Dai <- sqrt(diag(diag(AIS)))
  
  # image correlation matrix
  
  IR <- MASS::ginv(Dai)%*%IS%*%MASS::ginv(Dai)  
  
  AIR <- MASS::ginv(Dai)%*%AIS%*%MASS::ginv(Dai)  
  
  pcor <- AIR
  
  # negatives to partial correlation
  
  pcor[upper.tri(pcor)] <- -pcor[upper.tri(pcor)]
  
  pcor[lower.tri(pcor)] <- -pcor[lower.tri(pcor)]
  
  #return the partial correlation matrix
  
  return(pcor)
  
}
