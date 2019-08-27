#' @Title Function to prepare imputed data matrices for graphing with qgraph.
#' @name prep4graph
#' @author Alex daSilva
#' @return A list of lists containing properly formatted matrices, entries are 0 if the relationship is not significant
#' 
#' @param x a list of pooled imputations returned from 'mlVAR_rubin' or a list of analyzed data from 'mlVAR00'
#' @param type are the data pooled imputations (imputed) or a complete case analyzed by mlVAR00
#' 
#' @importFrom Matrix forceSymmetric

  
######need to add funtionality to get pcors and average; then update prep4graph!!!!!!!!

prep4graph <- function(x, type = c("imputed", "complete")) {
   
   graph_dataOut <- list()
   
   if (type == "imputed") {
     
     varNames <- names(x)[1:3]
     
   } else {
     
     varNames <- names(x[[1]])
     
   }
   
   #Grab the appropriate estimates (pcors for contemp, between; betas for temporal) to be used in graphing
   
   if (type == "imputed") {
      
      ests <- list(x[["temporal"]][["m_o_m"]], x[["between_pcors"]], x[["contemporaneous_pcors"]])
      
      ps <- list(x[["temporal"]][["pooled_p"]], x[["between-subjects"]][["pooled_p"]], x[["contemporaneous"]][["pooled_p"]])
      
   } else {
      
      ests <- list(x[[1]][["temporal"]][["means"]], x[[1]][["between-subjects"]][["pcor"]], x[[1]][["contemporaneous"]][["pcor"]])
      
      ps <- list(x[[1]][["temporal"]][["P-value"]], x[[1]][["between-subjects"]][["P-value"]], x[[1]][["contemporaneous"]][["P-value"]])
      
   }
   
  
   for (j in 1:length(varNames)) {
     
     # # get current network
     # 
     # if (type == "imputed") {
     #   
     #   ps <- x[[j]][["pooled_p"]]
     #   
     #   ests <- x[[j]][["m_o_m"]]
     #   
     # } else {
     #   
     #  ps <-  x[[1]][[j]][["P-value"]]
     #   
     #  ests <- x[[1]][[j]][["means"]] 
     #  
     # }
     
     # find sig relationships
     
     inds <- which(ps[[j]] < .05)
     
     out <- double()
     
     # if a relationship isnt sig, make it 0
     # if there are no sig relationships everything in the output vector 'out' is 0
     
     if (length(inds) != 0) {
       
       for (i in 1:length(c(ps[[j]]))) {
         
         if (i %in% inds) {
           
           out[i]  <- c(ests[[j]])[i]
           
         } else {
           
           out[i] <- 0
           
         }
         
       }
       
     } else {
       
       out <- rep(0, length(c(ps[[j]])))
       
     }
     
     #make a note if there are tiny beta coefs that are sig, if there are some this small make sure threshold to exclude is small than this in qgraph
     
     check <- (out > 0 & out < .00001) |  (out < 0 & out > -.00001)
     
     if (TRUE %in% c(check) == TRUE) {
       
       warning("check your minimum input into qgraph small beta values present")
       
     }
     
     #vector to matrix
     
     graphIN <- matrix(out, nrow = dim(ps[[j]])[1], dim(ps[[j]])[2], byrow = FALSE, dimnames = dimnames(ps[[j]])) 
     
     #force symmetric for qgraph, if contemp or btwn
     
     if (varNames[j] != "temporal") {
       
       graphIN <- Matrix::forceSymmetric(graphIN)
       
     } else {
       
       graphIN <- graphIN
       
     }
     
     graph_dataOut[[j]] <- t(graphIN)
     #graph_dataOut[[j]] <- graphIN
     
   }
   
   if (type == "imputed") {
     
     names(graph_dataOut) <- varNames
     
   } else {
     
     names(graph_dataOut) <- varNames
     
   }
   
   return(graph_dataOut)
   
 }