
win_btwn <- function(v, name){
  
  win <- scale(v, center = TRUE, scale = FALSE)
  
  btwn <- rep(mean(v, na.omit = TRUE), length(v))
  
  out <- cbind(win, btwn)
  
  colnames(out) <- c(paste("win", name, sep = "_"), paste("btwn", name, sep = "_") )
  
  return(out)
  
}
