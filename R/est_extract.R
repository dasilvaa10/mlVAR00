#' @Title Function to extract relevant coefficients from a lmer model.
#' @name est_extract
#' @author Alex daSilva
#' @return A dataframe of lmer model coefficients
#' 
#' @param mod a lme4 model


est_extract <- function(mod) {
  
  #grab coefs, ses, ts from model
  
  coefs <- as.data.frame(summary(mod)$coef)
  
  coefs$variable <- row.names(coefs)
  
  #move variable names from row.names to actual column
  
  coefs <- coefs[,c(4, 1, 2, 3)]
  
  row.names(coefs) <- NULL 
  
  return(coefs)
  
} 