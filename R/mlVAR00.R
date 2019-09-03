#' @Title Function to fit vector autoregressive models.
#' @name mlVAR00
#' @author Alex daSilva
#' @return A list of lists containing network results, formatted data, and underlying lme4 models
#' 
#' @param dat time-varying longituindal data to analyze
#' @param scale whether or not to standardized data
#' @param variables list of variables to analyze
#' @param ID name of id variable
#' @param temporal specifies random effects strucutre for first set of models. Can be "correlated", "interceptsOnly", or "orthogonal". Default is "correlated"
#' @param contemporaneous specifies random effects strucutre for second set of models. Can be "correlated", or "orthogonal". Default is "correlated"
#' @param timeVar name of temporal variable
#' @param timePoly option to add nth order polynomial terms for temporal variable
#' @param timeRandom should random effects be included for temporal variable
#' 
#' @importFrom MASS ginv
#' @importFrom lme4 lmer
#' @importFrom corpcor pseudoinverse


mlVAR00 <- function(dat, scale = FALSE, variables = NULL, ID = NULL, temporal = "correlated", contemporaneous = "correlated", timeVar = NULL, timePoly = NULL, timeRandom = FALSE, silenceLMER = TRUE ) {
  
  #variable specification and data formatting
  
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
  
  #if temporal variable is included, specify nth order polynomial effects
  
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
  
  formatted_data <- na.omit(formatted_data)
  
  #fit and first set of models; temporal and between network coefficients obtained
  
  temporalModels <- list()
  
  temporalModelsResids <- list()
  
  temporalSummary <- list()
  
  message("Check output for LMER convergence! Convergence issues may arise with a fully specified random effect stucture for a relatively larger number of variables")
  
  cat("\n", "estimating temporal and between-subject networks", "\n")
  
  pb <- txtProgressBar(min = 0, max = length(varNames), style = 3)
  
  for (i in 1:length(varNames)){
    
    #begin by specifying random effect structure
    
    trait2rm <- paste("btwn", varNames[i], sep = "_")
    
    #maximal random effect structure
    
    if (temporal == "correlated") {
      
      if (timeRandom == TRUE) {
        
        reffs <- colnames(formatted_data)[grep("win", colnames(formatted_data))]
        
        reffs <- c(reffs, timeRE)
        
      } else {
        
        reffs <- colnames(formatted_data)[grep("win", colnames(formatted_data))]
        
      }  
      
      reffsFrm <- paste("+(",paste(reffs, collapse = "+"), "|ID)", sep = "")
      
      #random intercept only
      
    } else if (temporal == "interceptOnly") {
      
      reffsFrm <- "+(1|ID)"
      
      #orthogonal random effects
      
    } else if (temporal == "orthogonal") {
      
      if (timeRandom == TRUE){
        
        reffs <- colnames(formatted_data)[grep("win", colnames(formatted_data))]
        
        reffs <- c(reffs, timeRE)
        
      }  else {
        
        reffs <- colnames(formatted_data)[grep("win", colnames(formatted_data))]
        
      }
      
      #create formula for random effects portion of model
      
      reffs_temp <- paste(paste("(0 +", reffs), "| ID)")
      
      reffs_temp <- paste(reffs_temp, collapse = "+")
      
      reffs_temp <- paste(reffs_temp, ")", sep = "")
      
      reffsFrm <- paste("+ ((1|ID) +", reffs_temp)
      
    } else {
      
      stop("temporal parameter not valid")
      
    }
    
    # obtain the variables to be used as fixed effects
    
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
    
    #fixed effects portion of formula
    
    feffsFrm <- paste(feffs, collapse = "+")
    
    #set i-th variable as the dependent variable
    
    dep <- paste(varNames[i], "~", sep = "")
    
    #final complete formula to pass to lmer
    
    frmla <- as.formula(paste(dep, feffsFrm, reffsFrm))
    
    #silence warnings
    
    if (silenceLMER == TRUE) {
      
      mod <- suppressWarnings(suppressMessages(lme4::lmer(frmla, data = formatted_data, REML =  FALSE)))
      
    } else {
      
      mod <- lme4::lmer(frmla, data = formatted_data, REML = FALSE)
      
    }
    
    #extract coefficients, ts, and ses
    
    temporalSummary[[i]]<- est_extract(mod)
    
    #save residuals for next set of models
    
    modResids <- residuals(mod)
    
    temporalModelsResids[[i]] <- modResids
    
    #save model itself
    
    temporalModels[[i]] <- mod
    
    setTxtProgressBar(pb, i)
    
  }
  
  close(pb)
  
  names(temporalModels) <- varNames
  
  #mark <- length(varNames) - 1
  
  #ind <- mark * (mark + 1)
  
  btwn_coefs <- list()
  
  btwn_ses <- list()
  
  win_coefs <- list()
  
  win_ses <- list()
  
  #get beta estimates and standard errors, separate between temporal and between-subject networks
  
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
  
  #get t-stats and pvalues for between-subjects network
  
  btwnMeanMat <- dat2mat(btwn_coefs, ind = length(varNames), temporal = FALSE, varNames = varNames)
  
  btwnSESMat <- dat2mat(btwn_ses, ind = length(varNames), temporal = FALSE, varNames = varNames)
  
  btwnTmat <- btwnMeanMat / btwnSESMat
  
  btwnPmat <- 2*(1-pnorm(abs(btwnTmat)))
  
  between_subjects <- list(btwnMeanMat, btwnSESMat, btwnTmat, btwnPmat)
  
  names(between_subjects) <- c("means", "SES", "T", "P")
  
  #get t-stats and pvalues for temporal network
  
  temporalMeanMat <- dat2mat(win_coefs, ind = length(varNames), temporal = TRUE, varNames = varNames)
  
  temporalSESMat <- dat2mat(win_ses, ind = length(varNames), temporal = TRUE, varNames = varNames)
  
  temporalTmat <- temporalMeanMat / temporalSESMat
  
  temporalPmat <- 2*(1-pnorm(abs(temporalTmat)))
  
  temporal <- list(temporalMeanMat, temporalSESMat, temporalTmat, temporalPmat)
  
  names(temporal) <- c("means", "SES", "T-value", "P-value")
  
  #get cor and pcor for between-subjects network
  #from https://github.com/SachaEpskamp/mlVAR
  
  mu_SD <- sapply(temporalModels, function(x) attr(lme4::VarCorr(x)[[1]], "stddev")[1])  #gets the intercept stddev
  
  D <- diag(1/mu_SD^2)
  
  inv <- D %*% (diag(length(varNames)) - btwnMeanMat)
  
  inv <- (inv + t(inv))/2
  
  inv <- forcePositive(inv)
  
  btwn_cov <- corpcor::pseudoinverse(inv)
  
  btwn_prec <- inv
  
  btwn_cor <- cov2cor(btwn_cov)
  
  btwn_pcor <- cor2pcor(btwn_cor, type = "cor")
  
  #store betas, ses, ts, ps, cors, and pcors in a list
  
  between_subjects <- list(btwnMeanMat, btwnSESMat, btwnTmat, btwnPmat, btwn_cor, btwn_pcor)
  
  names(between_subjects) <- c("means", "SES", "T-value", "P-value", "cor", "pcor")
  
  #estimate contemporaneous network
  #create data from residuals
  
  contempDat <- do.call("cbind", temporalModelsResids)
  
  colnames(contempDat) <- varNames
  
  contempDat <- data.frame(contempDat, ID = formatted_data$ID)
  
  contempModels <- list()
  
  contempSummary <- list()
  
  cat("\n","estimating contemporaneous network", "\n")
  
  pb <- txtProgressBar(min = 0, max = length(varNames), style = 3)
  
  #begin model fitting
  
  for ( i in 1:length(varNames)){
    
    #set variable i as dependent variable
    
    dep <- varNames[i]
    
    #get fixed effects
    
    fixeffs <- varNames[!(varNames %in% dep)]
    
    dep <- paste(dep, "~", sep = "")
    
    #formatting formulas for fixed effects
    #if only 2 variables used, no need to collapse
    
    if (length(fixeffs) <= 1) {
      
      fixeffsFrm <- fixeffs
      
    } else {
      
      fixeffsFrm <- paste(fixeffs, collapse = "+")
      
    }
    
    fixeffsFrm  <- paste("0", fixeffsFrm, sep  = "+")
    
    #maximal random effects
    
    if (contemporaneous == "correlated") {
      
      raneffsFrm <- paste(paste("(0", fixeffs, sep = "+"), "|ID)", sep = "")
      
      #orthogonal random effects
      
    } else if (contemporaneous == "orthogonal"){
      
      #formatting random effects formula
      
      reffs_temp <- paste(paste("(0 +", fixeffs), "| ID)")
      
      reffs_temp <- paste(reffs_temp, collapse = "+")
      
      reffs_temp <- paste(reffs_temp, ")", sep = "")
      
      raneffsFrm <- paste("(", reffs_temp, sep = "")
      
    }
    
    #frmla <- as.formula(paste(paste(dep, fixeffsFrm, sep = ""), raneffsFrm, sep = "+"))
    
    #complete formula
    
    frmla <- as.formula(paste(paste(dep, fixeffsFrm, sep = ""), paste(raneffsFrm, collapse = "+"), sep = "+"))
    
    if (silenceLMER == TRUE) {
      
      mod <- suppressWarnings(suppressMessages(lme4::lmer(frmla, data = contempDat, REML = FALSE)))
      
    } else {
      
      mod <- lme4::lmer(frmla, data = contempDat, REML = FALSE)
      
    }
    
    #store model
    
    contempModels[[i]] <- mod
    
    #extract betas, ses, ts
    
    contempSummary[[i]]<- est_extract(mod)
    
    setTxtProgressBar(pb, i)
    
  }
  
  close(pb)
  
  names(contempModels) <- varNames
  
  contemp_coefs <-list()
  
  contemp_ses <-list()
  
  #grab betas and ses
  
  for ( i in 1:length(varNames)) {
    
    contemp_coefs[[i]] <-  contempSummary[[i]]$Estimate
    
    contemp_ses[[i]] <- contempSummary[[i]]$`Std. Error`
    
  }
  
  #get t-stats and pvalues for temporal network
  
  contempMeanMat <- dat2mat(contemp_coefs, ind = length(varNames), temporal = FALSE, varNames = varNames)
  
  contempSESMat <- dat2mat(contemp_ses, ind = length(varNames), temporal = FALSE, varNames = varNames)
  
  contempTmat <- contempMeanMat / contempSESMat
  
  contempPmat <- 2*(1-pnorm(abs(contempTmat)))
  
  #get cors and pcors for contemporaneous network
  # from https://github.com/SachaEpskamp/mlVAR
  
  D <- diag(1/sapply(contempModels,sigma)^2)
  
  inv <- D %*% (diag(length(varNames)) - contempMeanMat)
  
  inv <- (inv + t(inv))/2
  
  inv <- forcePositive(inv)
  
  contemp_cov <- corpcor::pseudoinverse(inv)
  
  contemp_prec <- inv
  
  contemp_cor <- cov2cor(contemp_cov)
  
  contemp_pcor <- cor2pcor(contemp_cor, type = "cor")
  
  #store betas, ses, ts, ps, cors, and pcors in list
  
  contemporaneous <- list(contempMeanMat, contempSESMat, contempTmat, contempPmat, contemp_cor, contemp_pcor )
  
  names(contemporaneous) <- c("means", "SES", "T-value", "P-value", "cor", "pcor" )
  
  #store results from 3 network analyses
  
  results <-list(temporal, between_subjects, contemporaneous)
  
  names(results) <- c("temporal", "between-subjects", "contemporaneous")
  
  #only saving 2 types of models as temporal/between are estimated in one shot
  
  models <- list(temporalModels, contempModels)
  
  names(models) <- c("temporal", "contemporaneous")
  
  #save results, models and the formatted data frame
  
  list2return <- list( results, models, formatted_data, contempDat)
  
  names(list2return) <- c("results", "models", "data", "residuals")
  
  return(list2return)
  
}
