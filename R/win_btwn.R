#' @Title Function to separate a numeric vector for a given participant into within and between subject components.
#' @name win_btwn
#' @author Alex daSilva
#' @return A matrix of within centered values and a mean repeated between subject value
#' 
#' @param x data matrix
#' @param name character, variable name 


win_btwn <- function(v, name){
  
  win <- scale(v, center = TRUE, scale = FALSE)
  
  btwn <- rep(mean(v, na.omit = TRUE), length(v))
  
  out <- cbind(win, btwn)
  
  colnames(out) <- c(paste("win", name, sep = "_"), paste("btwn", name, sep = "_") )
  
  return(out)
  
}
