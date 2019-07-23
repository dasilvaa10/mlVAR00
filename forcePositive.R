# from https://github.com/SachaEpskamp/mlVAR

forcePositive <- function(x){
  x <- (x + t(x))/2
  
  if (any(eigen(x)$values < 0)){
    
    return(x - (diag(nrow(x)) * min(eigen(x)$values)-0.001)) 
    
  } else {
    
    return(x)
    
  }
  
}
