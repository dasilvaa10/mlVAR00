rubin_combine <- function(data_list, m){
  
  
  current_network_out <- list()
  
  for (j in 1:3){
    
    ind <- dim(data_list[[1]]$temporal$means)[1]
    
    l <- ind*ind
    
    crnt_network_label <- names(data_list[[1]][1:3])
    
    temp <- lapply(data_list, `[[`, crnt_network_label[j])
    
    temp_means <- lapply(temp, `[[`, "means")
    
    temp_ses <- lapply(temp, `[[`, "ses")
    
    temp_variances <- lapply(temp_ses, function(x) x^2)
    
    var_specific_variances <- list()
    
    var_specific_means <- list()
    
    ###remove hard-coding
    
    for (i in 1:l) { 
      
      var_specific_variances[[i]] <- sapply(temp_variances, function(x) x[i])
      
      var_specific_means[[i]] <- sapply(temp_means, function(x) x[i])
      
    }  
    
    within_variance <- sapply(var_specific_variances, mean)
    
    btwn_variance <- sapply(var_specific_means, var)
    
    mean_of_means <- sapply(var_specific_means, mean)
    
    rubins_variance <- double()
    
    lambda <- double()
    
    dfold <- double()
    
    t_val <- double()
    
    pooled_p <- double()
    
    for (i in 1:l) {
      
      rubins_variance[i] <- within_variance[i] + ((1+1/m)*btwn_variance[i])
      
      lambda[i] <- (btwn_variance[i] + (btwn_variance[i]) / m) / (rubins_variance[i])
      
      dfold[i] <- (m - 1) / lambda[i]^2
      
      t_val[i] <- mean_of_means[i]/sqrt(rubins_variance[i])
      
      pooled_p[i] <- 2*pt(abs(t_val[i]), dfold[i], lower=FALSE)
      
    }
    
    rn <- row.names(data_list[[1]]$mlVar00_object$results$temporal$means)
    
    cn <- colnames(data_list[[1]]$mlVar00_object$results$temporal$means)
    
    t_mat <- matrix(c(t_val),ind,ind,byrow = FALSE)
    
    colnames(t_mat) <- cn
    
    row.names(t_mat) <- rn
    
    p_mat <- matrix(c(pooled_p),ind,ind,byrow = FALSE)
    
    colnames(p_mat) <- cn
    
    row.names(p_mat) <- rn
    
    current_network_out[[j]] <- list(t_mat, p_mat)
    
    names(current_network_out[[j]]) <- c("t_value","pooled_p")
    
  }
  
  
  names(current_network_out) <- c("temporal_network","contemporaneous_network", "between-subjects_network")
  
  return(current_network_out)
  
}