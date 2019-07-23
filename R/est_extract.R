est_extract <- function(mod) {
  
  coefs <- as.data.frame(summary(mod)$coef)
  
  coefs$variable <- row.names(coefs)
  
  coefs <- coefs[,c(4, 1, 2, 3)]
  
  row.names(coefs) <- NULL 
  
  return(coefs)
  
} 