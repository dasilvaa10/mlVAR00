mlVAR00 <- function(dat, scale = FALSE, variables = NULL, ID = NULL, temporal = "correlated", contemporaneous = "correlated", timeVar = NULL, timePoly = NULL, timeRandom = FALSE, silenceLMER = TRUE ) {
  
  if (ID %in% colnames(dat) == FALSE){
    
    stop("make sure ID variable is specified correctly")
    
  }
  
  if (anyNA(match(variables, colnames(dat))) ==  TRUE) {
    
    stop("make sure variables are correctly specified")
    
  }
  
  if ("win_" %in% colnames(dat) | "btwn_" %in% colnames(dat) == TRUE ) {
    
    stop("Column names cannot contain the following text: 'btwn_', 'win_' ")
    
  }
  
  ID_ind <- grep(ID, colnames(dat))
  
  colnames(dat)[ID_ind] <- "ID"
  
  if (is.null(timeVar) == FALSE){
    
    dat <- dat[, c(variables, "ID", timeVar)]
    
  }  else {
    
    dat <- dat[, c(variables, "ID")]
    
  }
  
  formatted_data <- format_mlVAR(dat, scale = scale, timeVar = timeVar)
  
  if (is.null(timeVar) == FALSE){
    
    varNames <- colnames(dat)[!(colnames(dat) %in% c("ID", timeVar))]
    
  } else {
    
    varNames <- colnames(dat)[!(colnames(dat) %in% c("ID"))]
    
  }
  
  ###poly calc
  
  if (is.null(timePoly) == FALSE){
    
    if (is.null(timeVar) == FALSE & timePoly > 1) {
      
      polyDat <- poly(formatted_data[["time"]], timePoly, raw = T)
      
      polyDat <- polyDat[, -1, drop = FALSE]
      
      polyDat_names <- paste( "time" , colnames(polyDat), sep = "_")
      
      colnames(polyDat) <- polyDat_names
      
      polyDat <- scale(polyDat)
      
      formatted_data <- cbind(formatted_data, polyDat )
      
    }
    
    timeRE <- colnames(formatted_data)[grep("time", colnames(formatted_data))]
    
  }  
  
  temporalModels <- list()
  
  temporalModelsResids <- list()
  
  temporalSummary <- list()
  
  message("Check output for LMER convergence! Convergence issues may arise with a fully specified random effect stucture for a relatively larger number of variables")
  
  cat("\n", "estimating temporal and between-subject networks", "\n")
  
  pb <- txtProgressBar(min = 0, max = length(varNames), style = 3)
  
  for (i in 1:length(varNames)){
    
    trait2rm <- paste("btwn", varNames[i], sep = "_")
    
    if (temporal == "correlated") {
      
      if (timeRandom == TRUE) {
        
        reffs <- colnames(formatted_data)[grep("win", colnames(formatted_data))]
        
        reffs <- c(reffs, timeRE)
        
      } else {
        
        reffs <- colnames(formatted_data)[grep("win", colnames(formatted_data))]
        
      }  
      
      reffsFrm <- paste("+(",paste(reffs, collapse = "+"), "|ID)", sep = "")
      
    } else if (temporal == "interceptOnly") {
      
      reffsFrm <- "+(1|ID)"
      
    } else if (temporal == "orthogonal") {
      
      if (timeRandom == TRUE){
        
        reffs <- colnames(formatted_data)[grep("win", colnames(formatted_data))]
        
        reffs <- c(reffs, timeRE)
        
      }  else {
        
        reffs <- colnames(formatted_data)[grep("win", colnames(formatted_data))]
        
      }
      
      reffs_temp <- paste(paste("(0 +", reffs), "| ID)")
      
      reffs_temp <- paste(reffs_temp, collapse = "+")
      
      reffs_temp <- paste(reffs_temp, ")", sep = "")
      
      reffsFrm <- paste("+ ((1|ID) +", reffs_temp)
      
    } else {
      
      stop("temporal parameter not valid")
      
    }
    
    feffs <- colnames(formatted_data)[!(colnames(formatted_data) %in% c(trait2rm, varNames, "ID"))]
    
    #feffs <- sort(feffs) commented out on 7/21 trying to figure out matrix ordering
    
    # if (is.null(timeVar) == FALSE){
    #   
    #   if (timePoly == 1){
    #     
    #     polyAdds <- "time"
    #     
    #   } else {
    #     
    #     polyAdds <- c("time", polyDat_names)
    #     
    #   }
    #   
    #   feffs <- c(feffs, polyAdds)
    #   
    # }
    
    feffsFrm <- paste(feffs, collapse = "+")
    
    dep <- paste(varNames[i], "~", sep = "")
    
    frmla <- as.formula(paste(dep, feffsFrm, reffsFrm))
    
    if (silenceLMER == TRUE) {
      
      mod <- suppressWarnings(suppressMessages(lme4::lmer(frmla, data = formatted_data, REML =  FALSE)))
      
    } else {
      
      mod <- lme4::lmer(frmla, data = formatted_data, REML = FALSE)
      
    }
    
    temporalSummary[[i]]<- est_extract(mod)
    
    modResids <- residuals(mod)
    
    temporalModelsResids[[i]] <- modResids
    
    temporalModels[[i]] <- mod
    
    setTxtProgressBar(pb, i)
    
  }
  
  close(pb)
  
  names(temporalModels) <- varNames
  
  mark <- length(varNames) - 1
  
  ind <- mark * (mark + 1)
  
  btwn_coefs <- list()
  
  btwn_ses <- list()
  
  win_coefs <- list()
  
  win_ses <- list()
  
  for (i in 1:length(varNames)){
    
    temp <- temporalSummary[[i]][-1, ]
    
    btwn_inds <- grep("btwn_", temp$variable) 
    
    #win_inds <- c(1:nrow(temp))[!(1:nrow(temp) %in% btwn_inds)]
    
    win_inds <- grep("win_", temp$variable)
    
    btwn_coefs[[i]] <- temp$Estimate[btwn_inds]
    
    btwn_ses[[i]] <- temp$`Std. Error`[btwn_inds]
    
    win_coefs[[i]] <- temp$Estimate[win_inds]
    
    win_ses[[i]] <- temp$`Std. Error`[win_inds]
    
  }  
  
  ### get between estimates
  
  btwnMeanMat <- dat2mat(btwn_coefs, ind = length(varNames), temporal = FALSE, varNames = varNames)
  
  btwnSESMat <- dat2mat(btwn_ses, ind = length(varNames), temporal = FALSE, varNames = varNames)
  
  btwnTmat <- btwnMeanMat / btwnSESMat
  
  btwnPmat <- 2*(1-pnorm(abs(btwnTmat)))
  
  between_subjects <- list(btwnMeanMat, btwnSESMat, btwnTmat, btwnPmat)
  
  names(between_subjects) <- c("means", "SES", "T", "P")
  
  ### same for temporal
  
  temporalMeanMat <- dat2mat(win_coefs, ind = length(varNames), temporal = TRUE, varNames = varNames)
  
  temporalSESMat <- dat2mat(win_ses, ind = length(varNames), temporal = TRUE, varNames = varNames)
  
  temporalTmat <- temporalMeanMat / temporalSESMat
  
  temporalPmat <- 2*(1-pnorm(abs(temporalTmat)))
  
  temporal <- list(temporalMeanMat, temporalSESMat, temporalTmat, temporalPmat)
  
  names(temporal) <- c("means", "SES", "T-value", "P-value")
  
  ### get cor and pcor for between
  
  ### from https://github.com/SachaEpskamp/mlVAR
  
  mu_SD <- sapply(temporalModels, function(x) attr(lme4::VarCorr(x)[[1]], "stddev")[1])  #gets the intercept stddev
  
  D <- diag(1/mu_SD^2)
  
  inv <- D %*% (diag(length(varNames)) - btwnMeanMat)
  
  inv <- (inv + t(inv))/2
  
  inv <- forcePositive(inv)
  
  btwn_cov <- corpcor::pseudoinverse(inv)
  
  btwn_prec <- inv
  
  btwn_cor <- cov2cor(btwn_cov)
  
  ###
  
  btwn_pcor <- cor2pcor(btwn_cor)
  
  between_subjects <- list(btwnMeanMat, btwnSESMat, btwnTmat, btwnPmat, btwn_cor, btwn_pcor)
  
  names(between_subjects) <- c("means", "SES", "T-value", "P-value", "cor", "pcor")
  
  ###Contemporaneous estimates
  
  contempDat <- do.call("cbind", temporalModelsResids)
  
  colnames(contempDat) <- varNames
  
  contempDat <- data.frame(contempDat, ID = formatted_data$ID)
  
  contempModels <- list()
  
  contempSummary <- list()
  
  cat("\n","estimating contemporaneous network", "\n")
  
  pb <- txtProgressBar(min = 0, max = length(varNames), style = 3)
  
  for ( i in 1:length(varNames)){
    
    dep <- varNames[i]
    
    fixeffs <- varNames[!(varNames %in% dep)]
    
    dep <- paste(dep, "~", sep = "")
    
    if (length(fixeffs) <= 1) {
      
      fixeffsFrm <- fixeffs
      
    } else {
      
      fixeffsFrm <- paste(fixeffs, collapse = "+")
      
    }
    
    fixeffsFrm  <- paste("0", fixeffsFrm, sep  = "+")
    
    if (contemporaneous == "correlated") {
      
      raneffsFrm <- paste(paste("(0", fixeffs, sep = "+"), "|ID)", sep = "")
      
    } else if (contemporaneous == "orthogonal"){
      
      reffs_temp <- paste(paste("(0 +", fixeffs), "| ID)")
      
      reffs_temp <- paste(reffs_temp, collapse = "+")
      
      reffs_temp <- paste(reffs_temp, ")", sep = "")
      
      raneffsFrm <- paste("(", reffs_temp, sep = "")
      
    }
    
    #frmla <- as.formula(paste(paste(dep, fixeffsFrm, sep = ""), raneffsFrm, sep = "+"))
    
    frmla <- as.formula(paste(paste(dep, fixeffsFrm, sep = ""), paste(raneffsFrm, collapse = "+"), sep = "+"))
    
    if (silenceLMER == TRUE) {
      
      mod <- suppressWarnings(suppressMessages(lme4::lmer(frmla, data = contempDat, REML = FALSE)))
      
    } else {
      
      mod <- lme4::lmer(frmla, data = contempDat, REML = FALSE)
      
    }
    
    contempModels[[i]] <- mod
    
    contempSummary[[i]]<- est_extract(mod)
    
    setTxtProgressBar(pb, i)
    
  }
  
  close(pb)
  
  names(contempModels) <- varNames
  
  contemp_coefs <-list()
  
  contemp_ses <-list()
  
  for ( i in 1:length(varNames)) {
    
    contemp_coefs[[i]] <-  contempSummary[[i]]$Estimate
    
    contemp_ses[[i]] <- contempSummary[[i]]$`Std. Error`
    
  }
  
  contempMeanMat <- dat2mat(contemp_coefs, ind = length(varNames), temporal = FALSE, varNames = varNames)
  
  contempSESMat <- dat2mat(contemp_ses, ind = length(varNames), temporal = FALSE, varNames = varNames)
  
  contempTmat <- contempMeanMat / contempSESMat
  
  contempPmat <- 2*(1-pnorm(abs(contempTmat)))
  
  ### from https://github.com/SachaEpskamp/mlVAR
  D <- diag(1/sapply(contempModels,sigma)^2)
  
  inv <- D %*% (diag(length(varNames)) - contempMeanMat)
  
  inv <- (inv + t(inv))/2
  
  inv <- forcePositive(inv)
  
  contemp_cov <- corpcor::pseudoinverse(inv)
  
  contemp_prec <- inv
  
  contemp_cor <- cov2cor(contemp_cov)
  
  contemp_pcor <- cor2pcor(contemp_cor)
  
  ###
  
  contemporaneous <- list(contempMeanMat, contempSESMat, contempTmat, contempPmat, contemp_cor, contemp_pcor )
  
  names(contemporaneous) <- c("means", "SES", "T-value", "P-value", "cor", "pcor" )
  
  results <-list(temporal, between_subjects, contemporaneous)
  
  names(results) <- c("temporal", "between-subjects", "contemporaneous")
  
  models <- list(temporalModels, contempModels)
  
  names(models) <- c("temporal", "contemporaneous")
  
  list2return <- list( results, models, formatted_data)
  
  names(list2return) <- c("results", "models", "data")
  
  return(list2return)
  
}
