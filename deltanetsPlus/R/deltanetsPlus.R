#' DeltaNeTS+ function
#' 
#' The main function for DeltaNeTS+ for inferring gene perturbations for each sample. 
#' @import glmnet foreach doParallel Matrix
#' 
#' @param data a deltanetsPlus object or a numeric matrix/data.frame of log2FC data. For log2FC data, each row represents a gene and each column represents a sample.
#' @param slope The slope matrix from time-series log2FC data. The number of rows should be the same as the number of rows in the log2FC data. If a deltanetsPlus object is used or the log2FC data only include steady-state samples, no need to provide the \code{slope} input.
#' @param grn an nX2 matrix or data.frame of TF-gene interactions (optional). The 1st column is TFs and the 2nd column is the corresponding targets. The symbol in \code{grn} should be the same format of the rowanames of log2FC data. If grn is not provided, Lasso regression will be used to infer a DeltaNeTS+ model. Otherwise, ridge regression is used.
#' @param perturbation "group"(default) or "individual". "group" means predicting over all perturbation for each group. "individual" means predicting sample-by-sample perturbations. If \code{perturbation} is null, no need to provde \code{group}.
#' @param group A factor of group indication. The length of group should be equal to the number of samples. If \code{group} is null but \code{perturbation}="group", the group information will be obatined from the deltanetsPlus object provided.
#' @param cv.method "kfolds"(default) or "loo"(leave-one-out). If "loo" is selected, the input of \code{kfolds} will be ignored.
#' @param kfolds The number of folds for k-fold cross validation. Default==10.
#' @param lambda A sequence of lambda values for glmnet regression. We recommend the default value = 10^seq(-2,5,length.out=100).
#' @param cv.opt "lambda.min"(default) or "lambda.1se". For \code{perturbation}=="group", "lambda.1se" may perform better.
#' @param gset A subset of genes when you want to calculate perturbation for only those genes.
#' @param A.return A logical value. If it's ture, A matrix (network) will be returned. The default is FALSE, due to a large dimension of A in typical.
#' @param par A logical value. If TRUE, multipe cores are used for parallel computing. The default is FALSE (no parallel computation).
#' @param numCores The number of cores for parallel computing. The default is 4.
#' 
#' @return a list of two matrices
#' \item{P}{The matrix of gene perturbations. Each row corresponds to a gene following the same order as the one in the log2FC data, while each column corresponds to the levels of \code{group}.}
#' \item{A}{The n*n matrix of the inferred gene regulatory network. The (i,j)-th element of the matrix corresponds to the regulatory impact of gene j on gene i. The rwos and columns of the matrix correspond to genes following the same order as the one in the log2FC data.}
#' 
#' @export
deltanetsPlus <- function(data, slope=NULL, grn=NULL, perturbation=c("group","individual"), group=NULL, cv.method=c("kfolds","loo"), kfolds=10, lambda=NULL, cv.opt = c("lambda.min","lambda.1se"), gset=NULL, A.return=FALSE, par=FALSE, numCores=4){
  
  
  require(glmnet)
  require(Matrix)
  require(progress)
  
  
  if(class(data)=="deltanetsPlus"){
    
    lfc = data@lfc
    
    if(is.null(group)){
      group = data@experiment
    }
    
    if(!is.na(data@slope[1])){
      slope = data@slope
    }
    
    rm(data)
    
    
  }else if(is.matrix(data) | is.data.frame(data)){
    
    lfc = data
    rm(data)
    
  }else{
    
    stop("'data' should be either an objecti with 'deltanetsPlus' class or a matrix of log2FC data.")
    
  }

  
  m_lfc <- ncol(lfc)
  perturbation = perturbation[1]
  cv.method = cv.method[1]
  cv.opt = cv.opt[1]
  if(is.null(lambda)){
    lambda=10^seq(-2,5,length.out=100)
  }
  

  if(is.null(group) | perturbation=="individual"){
    group = factor(1:ncol(lfc))
  }else{
    if(length(group)!=ncol(lfc)){
      stop("The length of 'group' should be equal to the number of columns in 'lfc'.")
    }
    group = factor(group, levels=unique(group))
  }
  
  
  
  ## remove genes with no informative log2FC data (zeor expression or no variance across the samples)
  gn0 = nrow(lfc)
 
  if(!is.null(slope)){
    
    gi.tot = which(!apply(lfc==0,1,all) & apply(cbind(lfc,slope),1,function(x) length(unique(x))) > 2)
    lfc = lfc[gi.tot,]
    slope = slope[gi.tot,]
    gn = length(gi.tot)
    
    
    if(nrow(slope)!=nrow(lfc)){
      stope("The number of rows in 'slope' should be equal to the number of rows in 'lfc'.")}
    
    m_slope = ncol(slope)
    
    ## A factor for adjusting slope to the scale of log2FC
    lfc.norm = apply(lfc,1,function(x) norm(x,type="2"))
    slope.norm = apply(slope,1,function(x) norm(x,type="2"))
    
    
    rho = (lfc.norm/sqrt(m_lfc))/(slope.norm/sqrt(m_slope))
    # rho = rep(max(lfc.norm/sqrt(m_lfc))/max(slope.norm/sqrt(m_slope)),gn)
    # rho = (lfc.norm/sqrt(apply(lfc!=0,1,sum)))/(slope.norm/sqrt(apply(slope!=0,1,sum)))
    
    scaling.factor = function(gi){
      sc.fac = c(rep(1,m_lfc), rep(rho[gi], m_slope))
      sc.fac[is.na(sc.fac)]=1
      sc.fac[is.infinite(sc.fac)]=1
      return(sc.fac)
    }
    
    X0 = rbind(t(lfc),t(slope))
    
    Im = diag(m_lfc)
    Im = do.call(cbind, lapply(levels(group),function(group.i) apply(subset(Im,select=which(group==group.i)),1,sum)))
    Im = rbind(Im, matrix(0,nrow = m_slope, ncol = ncol(Im)))
    
    
  }else{
    m_slope=0
    gi.tot = which(!apply(lfc==0,1,all) & apply(lfc,1,function(x) length(unique(x))) > 2)
    lfc = lfc[gi.tot,]
    gn = length(gi.tot)
    
    scaling.factor = function(gi){1}
    
    X0 <- t(lfc)
    Im <- diag(m_lfc)
    Im <- do.call(cbind, lapply(levels(group), function(group.i) apply(subset(Im,select=which(group==group.i)),1,sum)))
  }
  
  
  
  glist = rownames(lfc)
  if(is.null(glist)){
    glist <- 1:gn
  }
  
  rm(list=c("lfc","slope"))
  
  message(paste("Removed ", gn0-gn," redundant gene(s) in DeltaNeTS+ analysis...\n" ,sep=""))
  

  
  if(!is.null(grn)){
    
    grn = as.matrix(grn)
    
    ## check grn format (tf,tg)
    if(ncol(grn)!=2){
      stop("'grn' should contain two columns of gene-gene interactions. 1st column:regulators, 2nd column:targets")
    }
    
    message("'grn' information is being applied.\n")
    
    ## Filtering interactions of genes only existing in the given expression data
    grn = unique(grn) # remove duplicated interactions
    
    grn = grn[is.element(grn[,1],glist) & is.element(grn[,2],glist),]
    if(nrow(grn)==0){
      stop("No gene in 'grn' was found in the expression data. Please check the format of gene names is the same between 'grn' and 'lfc'.")
    }
    
    ## Convert gene names in grn to indices for lfc
    grn_idx = cbind(match(grn[,1],glist),match(grn[,2],glist))
    grn_idx = grn_idx[which(grn_idx[,1]-grn_idx[,2]!=0),]## removed self-loops
    
    
    ## A set of downstream genes in grn
    dgi = sort(unique(grn_idx[,2]))
    if(!is.null(gset)){
      
      gset.i = match(gset,glist)
      dgi = intersect(dgi, gset.i)
      
      if(length(dgi)==0){
        stop("No regulator was foudn for the given gene set.")
      }
      
      message(paste("DeltsNeTS+ will predict the perturbation of the selected genes (",length(dgi)," genes) in 'grn'..\n", sep=""))
      
    }else{
      
      message(paste("DeltsNeTS+ will predict the perturbation of total regulon (",length(dgi)," genes) in 'grn'..\n", sep=""))
      
    }

    
    grn_usage = TRUE
    alpha = 0 ## ridge regression
    
    
    
  }else{
    
    grn_usage = FALSE
    dgi = 1:nrow(lfc)
    
    if(!is.null(gset)){
      
      gset.i = match(gset,glist)
      dgi = intersect(dgi, gset.i)
      
      if(length(dgi)==0){
        stop("No regulator was foudn for the given gene set.")
      }
      
      message(paste("DeltsNeTS+ will predict the perturbation of the selected genes (",length(dgi)," genes) in 'grn'..\n", sep=""))
      
    }
    
    alpha = 1 ## lasso regression
    message("DeltaNeTS+ will use Lasso regression for computing gene perturbations..\n")
  }
 
  #--------------------------------------------------------------------------------
  # DeltaNeTS+ model learning
  #--------------------------------------------------------------------------------
  ## Progress bar
  pb <- progress_bar$new(total=length(dgi))
  
  if(par){
    message("Parallel computing is on..\n")
    
    
    ## DeltaNeTS+ analysis with multiple work clusters
    require(doParallel)
    require(foreach)

    
    cl <- makeCluster(numCores,outfile='')
    registerDoParallel(cl)
    
    
    
    Sys.sleep(3)
    Beta <- foreach(j=1:length(dgi),.combine = cbind,.packages = "glmnet")%dopar%{
      gene.j = glist[dgi[j]]
      pb$tick(1)
      Sys.sleep(1 / 100)}
      
    
      if(grn_usage){
        tfi = grn_idx[which(grn_idx[,2]==dgi[j]),1]
        
      }else{
        tfi <- setdiff(1:gn, dgi[j]) ## all genes except itsef as regulators
      }
      
      X =subset(X0,select=tfi)*scaling.factor(dgi[j])
      x.norm = apply(X, 2, function(xi) norm(xi,type="2"))
      x.norm[x.norm==0]=1
      X = t(t(X)/x.norm)
      X = cbind(X, Im)
      
      
      y <- X0[,dgi[j]]*scaling.factor(dgi[j])
      y.nzi = which(y!=0)
      y = y[y.nzi]
      X = X[y.nzi,]
      x.nzi = which(apply(X!=0,2,sum)>0)
      
      if(length(unique(y)) < 3){
        r <- Matrix(0, ncol=1, nrow=(gn+ncol(Im)), sparse=TRUE)
      }else{
        
        if(cv.method=="loo"){
          kfolds = length(y.nzi)
        }
        ## model prediction using glmnet
        cvfit <- cv.glmnet(X[,x.nzi], y, family='gaussian', alpha = alpha, 
                           intercept = FALSE, nfolds = kfolds, standardize=FALSE,
                           lambda=lambda)
        
        # coef.cvfit <- (coef(cvfit, s='lambda.1se')+coef(cvfit, s='lambda.min'))/2
        coef.cvfit <- coef(cvfit, s=cv.opt)
        # bj <- rep(0,ncol(X))
        bj <- rep(0, length(x.nzi))
        bj[coef.cvfit@i] <- coef.cvfit@x
        r <- Matrix(0, ncol=1, nrow=(gn+ncol(Im)), sparse=TRUE)
        r[c(tfi,(gn+1):(gn+ncol(Im)))[x.nzi]] <- bj

      }
      
      r
      
 
    }else{
    ## No parallel computing
    Sys.sleep(3)
    Beta <- Matrix(0,nrow = (gn+ncol(Im)), ncol = length(dgi), sparse=TRUE)
    
    for (j in 1:length(dgi)){
      pb$tick(1)
      Sys.sleep(1 /100)

      if(grn_usage){
        tfi = grn_idx[which(grn_idx[,2]==dgi[j]),1]
        
      }else{
        tfi <- setdiff(1:gn, dgi[j]) ## all genes except itsef as regulators
      }
      
      X = subset(X0,select=tfi)*scaling.factor(dgi[j])
      x.norm = apply(X, 2, function(xi) norm(xi,type="2"))
      x.norm[x.norm==0]=1
      X = t(t(X)/x.norm)
      X = cbind(X, Im)
      
      
      y <- X0[,dgi[j]]*scaling.factor(dgi[j])
      y.nzi = which(y!=0)
      y = y[y.nzi]
      X = X[y.nzi,]
      x.nzi = which(apply(X!=0,2,sum)>0)
      
      if(length(unique(y)) < 3){
        Beta[,j] <- 0
      }else{
        
        if(cv.method=="loo"){
          kfolds = length(y.nzi)
        }
        
        cvfit <- cv.glmnet(X[,x.nzi], y, family='gaussian', alpha = alpha, 
                           intercept = FALSE, nfolds = kfolds, standardize=FALSE,
                           lambda=lambda)
        
        # coef.cvfit <- (coef(cvfit, s='lambda.1se')+coef(cvfit, s='lambda.min'))/2
        coef.cvfit <- coef(cvfit, s=cv.opt)
        # bj <- rep(0,ncol(X))
        bj <- rep(0, length(x.nzi))
        bj[coef.cvfit@i] <- coef.cvfit@x
        Beta[c(tfi,(gn+1):(gn+ncol(Im)))[x.nzi],j] <- bj
      }

    }
  }

  message("\n")
  message("Returning DeltaNeTS+ prediction..\n")
  if(A.return){
    A <- Matrix(0, nrow=gn, ncol=gn, sparse=TRUE)
    A[dgi,] <- Matrix::t(Beta[1:gn,])
    rownames(A) = glist
    colnames(A) = glist
    }else{A=NULL}
    
  
    P <- Matrix(0, nrow=gn, ncol=ncol(Im), sparse=TRUE)
    P[dgi,] <- Matrix::t(Beta[(gn+1):(gn+ncol(Im)),])
    colnames(P) = levels(group)
    rownames(P) = glist

  return(list(P=P,A=A))

}