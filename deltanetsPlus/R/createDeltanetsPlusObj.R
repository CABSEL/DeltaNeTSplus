#' creating a deltanetsPlus object.
#'
#'  This function creastes a deltanetsPlus object, which is an input format of \code{deltanetsPlus} function.
#' @param lfc The numeric matrix or data.frame of log2FC data. Each row represents a gene and each column represents a sample.
#' @param pval p-values for log2FC values. If pval==NULL, no filtering of log2FC data.
#' @param tp A vector of time points. The length of the vector should be the same as the number of samples. If steady-state samples, set the corresponding elements of tp to zero.
#' @param experiment A numeric vector for indicating the sampe experiment/perturbation. The (time-series) samples from the same experiment/perturbation should have the same unique index. The length of the vector should be the same as the number of samples.
#' @param threshold A threshold value for \code{pval}.
#'
#' @return a deltanetsPlus object
#'
#' @export
createDeltanetsPlusObj <- function(lfc=lfc, pval=NULL, tp=NULL, experiment=NULL, p.thres=0.05){
  require(methods)

  if(is.null(tp)){

    oi = 1:ncol(lfc)
    meta = data.frame(tp=rep(0,ncol(lfc)), experiment=1:ncol(lfc))
    experiment.type = rep("ss",ncol(lfc))
    message("No time information was provided.\nAll samples will be considered steady-state samples.")

  }else if(!is.numeric(tp)){

    stop("'tp' should be a numeric vector with the same length of ncol(lfc).")

  }else if(length(tp)!=ncol(lfc) | length(experiment)!=ncol(lfc)){

    stop("The length of 'tp' and 'experiment' should be equal to ncol(lfc).")

  }else if(!all(dim(cbind(tp,experiment))==dim(unique(cbind(tp,experiment))))){

    stop("Duplicated time points in the same experimental group. Each time point should be unique in an experiment.")

  }else{

    meta = data.frame(tp=tp, experiment=experiment)
    oi = order(meta$experiment, meta$tp)
    meta = meta[oi,]
    experiment.type = unlist(lapply(unique(meta$experiment),function(expi){nExp=sum(meta$experiment==expi); if(nExp>1){rep("ts",nExp)}else{"ss"}}))
  }




  if(!is.null(pval)){

    if(!all(dim(lfc)==dim(pval))){
      stop("For lfc filtering, the dimension of 'pval' should be equal to the dimension of 'lfc'.")
    }


    if(p.thres>1 | p.thres<0){
      stop("The p-value threshold 'p.thres' should be between 0 and 1 (default=0.05).")
    }


    ##  Filter and interpolate gene expressions with statistical significance
    lfc.new = unsig2interp(lfc[,oi], pval[,oi], meta$experiment, meta$tp, threshold=p.thres)
    filt=TRUE

  }else{
    lfc.new = lfc[,oi]
    filt = FALSE
  }




  if(sum(experiment.type=="ts")>0){

    slope <- generateSlope(lfc = lfc.new, tp = meta$tp, group = meta$experiment)
    rownames(slope) = rownames(lfc.new)

  }else{
    slope=NULL
  }

  ## Back to original sample order
  oi2 = order(oi)
  lfc.new=lfc.new[,oi2]
  experiment.type=experiment.type[oi2]
  experiment = experiment[oi2]

  ## create deltanetsPlus object
  setClass("deltanetsPlus", representation(lfc = "matrix", slope = "matrix",
                                       experiment = "factor",
                                       experiment.type="character",
                                       filt="logical"),
           prototype(slope = matrix()),
           where= topenv())

  setMethod("show","deltanetsPlus",function(object){
    cat("size lfc: ")
    print(dim(object@lfc));
    cat("size slope: ")
    if(!is.na(object@slope)[1]){
      print(dim(object@slope))
    }else{print(NULL)}

    cat("lfc filtering/interpolation:",object@filt,"\n")
    cat("Total number of experiments:", length(levels(object@experiment)),"\n")
    cat("Number of time-series samples:", sum(object@experiment.type=="ts"),"\n")
    cat("Number of steady-state samples:", sum(object@experiment.type=="ss"),"\n")
  })



  setMethod(f="cbind2",signature=c(x = "deltanetsPlus", y = "deltanetsPlus"),definition=function(x, y){
              lfc.comb = cbind2(x@lfc,y@lfc)
              if(is.na(x@slope[1])){x.slope=NULL}else{x.slope=x@slope}
              if(is.na(y@slope[1])){y.slope=NULL}else{y.slope=y@slope}
              slope.comb = cbind2(x.slope, y.slope)
              experiment.comb = c(as.character(x@experiment), as.character(y@experiment))
              experiment.comb = factor(experiment.comb,levels=unique(experiment.comb))
              experiment.type.comb = c(x@experiment.type, y@experiment.type)
              filt.comb = (x@filt | y@filt)

              obj.comb = new("deltanetsPlus", lfc = lfc.comb, slope=slope.comb, experiment=factor(experiment.comb,levels=unique(experiment.comb)), experiment.type=experiment.type.comb, filt=filt.comb)
              return(obj.comb)
            })



  if(!is.null(slope)){
    d <- new("deltanetsPlus", lfc = lfc.new, slope=slope, experiment=factor(experiment,levels=unique(experiment)), experiment.type=experiment.type,filt=filt)
  }else{
    d <- new("deltanetsPlus", lfc = lfc.new, experiment=factor(experiment,levels=unique(experiment)),experiment.type=experiment.type,filt=filt)
  }



  if(filt){
    message("lfc was filtered/interpolated using 'pval'.")
  }
  if(!is.na(d@slope[1])){
    message("Slope data were included.")
  }

  return(d)
}
