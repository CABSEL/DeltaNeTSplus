#' rank genes from P matrix
#' 
#' The function for ranking genes per sample, based on absolute magnitudes in estimates of perturbation impacts
#' @param scoremat a gene score matrix returned by \code{deltanetsPlus} function.
#' @param glist The list of genes corresponding to the rows of P matrix. If glist is not given, gene symbols are just numbers (e.g. G1, G2, G3, ...)
#' @return a list of rankOfGenes and rankedGenes
#' \item{rankOfGenes}{The matrix of gene ranks. Rows and columns correspond to genes and samples in the same order as the one in P matrix, respectively.A maximum rank is assigned for the genes with zero values in the P matrix}
#' \item{rankedGenes}{A list of ranked gene lists. Each element in the list includes a ranked gene list for the corresponding sample. The genes closer to the front are the gene targets with the higher confidence. }
#' 
#' @export
rankp <- function(scoremat, glist=NULL){

  n = nrow(scoremat)
  
  if(is.null(glist)){glist = paste("G",1:n,sep="") }

rankOfGenes = apply(abs(scoremat),2,function(x) order(order(x,decreasing=TRUE)))

rankedGenes = apply(abs(scoremat),2, function(x) glist[order(x,decreasing=TRUE)])

return(list(rankOfGenes=rankOfGenes, rankedGenes=rankedGenes))
}

