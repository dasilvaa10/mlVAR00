# Function in part by by G. Jay Kerns, Ph.D., Youngstown State University (http://tolstoy.newcastle.edu.au/R/e2/help/07/08/22816.html)

cor2pcor <- function(cormat) {
  
  X <- cormat 
  
  iX <- MASS::ginv(X) 
  
  S2 <- diag(diag((iX^-1)))
  
  AIS <- S2%*%iX%*%S2 # anti-image covariance matrix
  
  IS <- X+AIS-2*S2 # image covariance matrix
  
  Dai <- sqrt(diag(diag(AIS)))
  
  IR <- MASS::ginv(Dai)%*%IS%*%MASS::ginv(Dai)  # image correlation matrix
  
  AIR <- MASS::ginv(Dai)%*%AIS%*%MASS::ginv(Dai)  
  
  pcor <- AIR
  
  pcor[upper.tri(pcor)] <- -pcor[upper.tri(pcor)]
  
  pcor[lower.tri(pcor)] <- -pcor[lower.tri(pcor)]
  
  return(pcor)
  
}
