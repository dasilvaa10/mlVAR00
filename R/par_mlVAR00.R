#' @Title Function to fit a mlVAR00 model and extract coefficients to later pool under multiple imputation framework .
#' @name par_mlVAR00
#' @author Alex daSilva
#' @return A list of lists containing network coefficients, ses, and a mlVAR00 object
#' 
#' @param dat time-varying longituindal data to analyze
#' @param variables list of variables to analyze
#' @param ID name of id variable
#' @param rfStructure list specifing random effects strucutre for both sets of models. Can be "correlated", "interceptsOnly", or "orthogonal". Default is "correlated"
#' @param timeArgs list specifying temporal related arguments i.e., list(timeVar, timePoly, timeRandom)
#' @param scale whether or not to standardized data

#' @importFrom MASS ginv
#' @importFrom lme4 lmer
#' @importFrom corpcor pseudoinverse

par_mlVAR00 <- function(dat, variables, ID, rfStructure = list("correlated", "correlated"), timeArgs = list(NULL, NULL, FALSE), scale = TRUE ){
  
  #fit model
  
  fit <- mlVAR00(dat, variables = variables, ID = ID, temporal = rfStructure[[1]], contemporaneous = rfStructure[[2]], scale = scale, timeVar = timeArgs[[1]], timePoly = timeArgs[[2]], timeRandom = timeArgs[[3]])
  
  # extract temporal estimates
  
  temporal_means <- fit$results$temporal$means
  
  temporal_se <- fit$results$temporal$SES
  
  temporal_dat <- list(temporal_means, temporal_se)
  
  names(temporal_dat) <- c("means","ses")
  
  #extract contemporaneous estimates
  
  contemporaneous_means <- fit$results$contemporaneous$means
  
  contemporaneous_se <- fit$results$contemporaneous$SES
  
  contemporaneous_dat <- list(contemporaneous_means, contemporaneous_se)
  
  names(contemporaneous_dat) <- c("means","ses")
  
  #extract between-subject estimates
  
  between_means <- fit$results$`between-subjects`$means
  
  between_se <- fit$results$`between-subjects`$SES
  
  between_dat <- list(between_means, between_se)
  
  names(between_dat) <- c("means","ses")
  
  networks <- list(temporal_dat, contemporaneous_dat, between_dat, fit)
  
  names(networks) <- c("temporal", "contemporaneous", "between-subjects", "mlVar00_object")
  
  return(networks)
  
}

