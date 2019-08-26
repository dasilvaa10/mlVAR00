#' @Title Function to prepare imputed data matrices for graphing with qgraph.
#' @name prep4graph
#' @author Alex daSilva
#' @return A list of lists containing properly formatted matrices, entries are 0 if the relationship is not significant
#' 
#' @param x a list of pooled imputations returned from 'mlVAR_rubin'
#' @importFrom Matrix forceSymmetric


prep4graph < - function(x) {
   
   graph_dataOut <- list()
   
   for (j in 1:length(x)) {
     
     # get current network
     
     ps <- x[[j]][["pooled_p"]]
     
     # find sig relationships
     
     inds <- which(ps < .05)
     
     ests <- x[[j]][["m_o_m"]]
     
     out <- double()
     
     # if a relationship isnt sig, make it 0
     # if there are no sig relationships everything in the output vector 'out' is 0
     
     if (length(inds) != 0) {
       
       for (i in 1:length(c(ps))) {
         
         if (i %in% inds) {
           
           out[i]  <- c(ests)[i]
           
         } else {
           
           out[i] <- 0
           
         }
         
       }
       
     } else {
       
       out <- rep(0, length(c(ps)))
       
     }
     
     #make a note if there are tiny beta coefs that are sig, if there are some this small make sure threshold to exclude is small than this in qgraph
     
     check <- (out > 0 & out < .00001) |  (out < 0 & out > -.00001)
     
     if (TRUE %in% c(check) == TRUE) {
       
       warning("check your minimum input into qgraph small beta values present")
       
     }
     
     #vector to matrix
     
     graphIN <- matrix(out, nrow = dim(ps)[1], dim(ps)[2], byrow = TRUE, dimnames = dimnames(ps)) 
     
     #force symmetric for qgraph
     
     graphIN <- Matrix::forceSymmetric(graphIN)
     
     graph_dataOut[[j]] <- t(graphIN)
     
   }
   
   names(graph_dataOut) <- names(x)
   
   return(graph_dataOut)
   
 }