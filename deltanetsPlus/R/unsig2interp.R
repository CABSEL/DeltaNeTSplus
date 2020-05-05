#' interpolating log2FC values with p-value > threshold
#' 
#' The function for choosing log2FC values with p-value < threshold and interpolating the values for the time points between them.
#' @param lfc log2FC data
#' @param lfc.pval p-values for log2FC values, e.g. which can be obtained when averaging replicates using limma.
#' @param tp A vector of time points of the samples in the matrix lfc. The length of the vector should be the same as the number of samples (i.e. the number of columns in the matrix \code{lfc}). 
#' @param group A vector of indices indicating the set of samples from a particular drug/compound treatment. The (time-series) samples from the same drug treatment experiment should have the same unique index. The length of the vector should be the same as the number of samples.
#' @param threshold A threshold value for \code{lfc.pval}.
#' @return \item{lfc.new}
#' @export
unsig2interp <- function(lfc, lfc.pval, group, tp, threshold){
  
  lfc.trim = lfc
  lfc.trim[lfc.pval > threshold] = 0
  
  lfc.new = numeric(0)
  for(k in unique(group)){
    gri = which(group==k)
    if(length(gri)>1){
      
        lfc_k = lfc.trim[,gri]
        tp_k = tp[gri]
        
        lfc_k.last = apply(lfc_k, 1, function(x){nz.ind=which(x!=0); if(length(nz.ind)>0){x[nz.ind[which.max(nz.ind)]]}else{0}})
        
        lfc_k[,ncol(lfc_k)] = lfc_k.last
        
        lfc_k.interp = do.call(rbind, lapply(1:nrow(lfc_k), function(i){ind = unique(c(1, which(lfc_k[i,]!=0)));if(length(ind)>1){yout = approx(x=tp_k[ind], y=lfc_k[i,ind], xout=tp_k)$y; return(yout)}else{rep(0,ncol(lfc_k))}}))
    }else{
      lfc_k.interp = lfc.trim[,gri]
    }
    
    lfc.new = cbind(lfc.new, lfc_k.interp)
  }

  rownames(lfc.new) = rownames(lfc)
  colnames(lfc.new) = colnames(lfc)
  return(lfc.new)  
}
