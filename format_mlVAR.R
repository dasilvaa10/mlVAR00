format_mlVAR <- function(dat, scale = FALSE, timeVar = NULL){
  
  if (is.null(timeVar) == FALSE){
    
    time_ind <- grep(timeVar, colnames(dat))
    
    times<- dat[,c("ID", timeVar )]
    
    dat <- dat[, -c(time_ind)]
    
    timesSplit <- split.data.frame(times, times$ID)
    
    for (i in 1:length(timesSplit)) {
      
      lastRow <- nrow(timesSplit[[i]])
      
      timesSplit[[i]] <- timesSplit[[i]][-lastRow, ]
      
    }
    
  }
  
  if (scale == TRUE) {
    
    id_ind <- grep("ID", colnames(dat))
    
    dat[, - id_ind] <- scale(dat[, - id_ind])
    
  }
  
  varNames <- colnames(dat)[!(colnames(dat) %in% c("ID"))]
  
  dat_split <- split.data.frame(dat, dat$ID)
  
  subDat <- list()
  
  for (i in 1:length(dat_split)){
    
    temp <-list()
    
    for (j in 1:length(varNames)){
      
      temp[[j]] <- win_btwn(dat_split[[i]][,j], varNames[j])
      
    }
    
    dat_temp <- do.call("cbind", temp)
    
    dat_temp2 <- dat_temp[1:(nrow(dat_temp)- 1), ]
    
    if (is.null(timeVar) == TRUE) {
      
      subDat[[i]] <- cbind(dat_split[[i]][-1, c(varNames)], dat_temp2, ID = rep(unique(dat_split[[i]]$ID), nrow(dat_temp2)))
      
    } else {
      
      subDat[[i]] <- cbind(dat_split[[i]][-1, c(varNames)], dat_temp2, ID = rep(unique(dat_split[[i]]$ID), nrow(dat_temp2)), time = timesSplit[[i]][, 2] )
      
    }  
    
  }
  
  datReformat <- do.call("rbind", subDat)
  
  return(datReformat)
  
}
