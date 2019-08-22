#' @Title Function to randomly create missing data among multilevel data.
#' @name ml_missing
#' @author Alex daSilva
#' @return A dataset containing missing values
#' 
#' @param dat data set to create missingness in
#' @param variables variables to add missingness to
#' @param temporalVar time variable
#' @param ID id variable
#' @param ranSeq a sequence of specifying a range of values to create missinging from by ID, by variable


ml_missing <- function(dat, variables, temporalVar, ID, ranSeq) { 
  
  dat <- dat[, c(ID, variables, temporalVar)]
  
  # split data into list of lists, 1 for each ID
  
  datSplit <- split.data.frame(dat, dat[[ID]])
  
  nPerson <- length(datSplit)
  
  iList <- list()
  
  for (i in 1:nPerson){
    
    jList <- list()
    
    nTime <- unique(sapply(datSplit[[i]], nrow))
    
    # for variable (on the i-th subject) randomly remove a portion of the data vector
    
    for (j in 1:length(variables)) {
      
      perMiss <- sample(ranSeq, 1)
      
      nMiss<- perMiss * nTime
      
      missInds <- sort(sample(1:nTime, nMiss, replace = FALSE ))
      
      tempVec <- datSplit[[i]][, variables[j] ] 
      
      tempVec[missInds]<- NA
      
      jList[[j]] <-  tempVec
      
    }
    
    # combine the columns containing missing data
    
    tempDat <- do.call("cbind", jList)
    
    colnames(tempDat) <- variables
    
    # add back the id and temporal variable
    
    tempDat <- data.frame(tempDat, ID = datSplit[[i]][[ID]] )
    
    tempDat[[temporalVar]] <- datSplit[[i]][[temporalVar]]
    
    iList[[i]]<- tempDat
    
  }
  
  dat_cmbnd <- do.call("rbind", iList)
  
  return(dat_cmbnd)
  
}