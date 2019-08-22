#' @Title Function to pool mlVAR models under multiple imputation framework .
#' @name mlVAR_rubin
#' @author Alex daSilva
#' @return A list of lists containing pooled network data
#' 
#' @param data_list list of coefficients from network models
#' @param m number of imputed data sets, equal to the length of data_list


mlVAR_rubin <- function(data_list, m){
  
  current_network_list <- list()
  
  for (j in 1:3){
    
    #extract first network, here temporal
    
    current_NW <- names(data_list[[1]])[1:3]
    
    temp_nets <- lapply(data_list , `[[` , current_NW[j])
    
    #get first dim & row/ column names
    
    ind <- dim(temp_nets[[1]]$means)[1]
    
    ind2 <- ind^2
    
    rn <- row.names(temp_nets[[1]]$means)
    
    cn <- colnames(temp_nets[[1]]$means)
    
    #get a list of the means and ses
    
    means <- lapply(temp_nets, `[[`, "means" )
    
    ses <- lapply(temp_nets, `[[`, "ses" )
    
    # get 40 means and ses for the first variable e.g [1,1] or [1] and continue to pooled p_value
    # return beta value, tvals, rubin se, and pooled p values
    
    m_o_m <- double()
    
    t_vals <- double()
    
    rub_se <- double()
    
    pooled_p <- double()
    
    for (i in 1:ind2){
      
      current_mean <- sapply(means, `[[`, i)
      
      current_var <- sapply(ses, `[[`, i)^2
      
      # between imputation variance
      
      btwn_variance <- var(current_mean)
      
      # withing imputation variance
      
      win_variance <- mean(current_var)
      
      # pooled variance
      
      total_var <- win_variance + btwn_variance + (btwn_variance/m)
      
      # pooled point estimate
      
      meanOFmeans <- mean(current_mean)
      
      t_pooled <- meanOFmeans / sqrt(total_var)
      
      #calculate degrees of freedom
      
      lambda  <- (btwn_variance + (btwn_variance / m)) / total_var
      
      dfOld <- (m - 1) / lambda^2
      
      pooled_pval <- 2*pt(abs(t_pooled), dfOld, lower=FALSE)
      
      m_o_m[i] <-  meanOFmeans
      
      t_vals[i] <- t_pooled
      
      rub_se[i] <- sqrt(total_var)
      
      pooled_p[i] <- pooled_pval
      
    }
    
    # matrices of information to return
    
    m_o_m_mat <- matrix(m_o_m, nrow = ind, ncol = ind, dimnames = list(rn,cn))
    
    t_mat <- matrix(t_vals, nrow = ind, ncol = ind, dimnames = list(rn,cn))
    
    rub_mat <- matrix(rub_se, nrow = ind, ncol = ind, dimnames = list(rn,cn))
    
    p_mat <- matrix(pooled_p, nrow = ind, ncol = ind, dimnames = list(rn,cn))
    
    current_network_list[[j]] <- list(m_o_m_mat, t_mat, rub_mat, p_mat)
    
    names(current_network_list[[j]]) <- c("m_o_m","t_value","rub_se", "pooled_p")
    
  }
  
  names(current_network_list) <- c("temporal", "contemporaneous", "between")
  
  return(current_network_list)
  
}
