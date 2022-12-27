#
#A Fast and Efficient Approach for Gene-based Association Studies of Ordinal Phenotypes
#
#OR_ACAT_CR:association test between rare and common variants and an ordinal phenotype 
#
#ord_data_pv:variant-level p values 
#
#order_Burden:Burden method to aggregate rare variants
#
#ACAT:Aggregated Cauchy Association Test 
#a general, powerful, and computationally efficient P-value combination method 
#(https://github.com/yaowuliu/ACAT)
#
#x1 :genotype score
#ord_y:ordinal phenotype
#xc:covariate
#rvt:rare variant threshold
#################################
OR_ACAT_CR<-function(G=x1,obj=ord_y,xc,Cweights.beta=c(1,25),Rweights.beta=c(1,25),common.weights=NULL,rare.weights=NULL,rvt=rvt){
  G=as.matrix(G)
  obj=as.matrix(obj)
  xc=as.matrix(xc)
  
  ### extract rare variant
  maf = colMeans(G)/2 
  rare.pre=which(maf<=rvt & maf>0)
  rare.G=G[,rare.pre]
  
  ### extract common variant
  common.pre=which(maf>rvt & maf<0.5)
  common.G=G[,common.pre]
  
  ### check weights
  if (!is.null(common.weights) && length(common.weights)!=ncol(common.G)){
    stop("The length of common weights must equal to the number of common variants!")
  }
  if (!is.null(rare.weights) && length(rare.weights)!=ncol(rare.G)){
    stop("The length of rare weights must equal to the number of rare variants!")
  }
  
  
  ###check G
  if (length(maf)==0){
    stop("The genotype matrix do not have non-zero element!")
  }

  ### rare-variant p-values and weights
  if (!is.null(length(rare.G))){
    
    if (is.null(ncol(rare.G))){
      if (is.null(rare.weights)){
        rare.mafs=mean(rare.G)
        rare.W<-(dbeta(rare.mafs,Rweights.beta[1],Rweights.beta[2])*sqrt(rare.mafs*(1-rare.mafs)))^2
      }else{
        rare.W<-mean(rare.weights)
      }
    }else{
      if (is.null(rare.weights)){
        rare.MAF=colMeans(rare.G)/2
        rare.mafs=mean(rare.MAF)
        rare.W<-(dbeta(rare.mafs,Rweights.beta[1],Rweights.beta[2])*sqrt(rare.mafs*(1-rare.mafs)))^2
      }else{
        rare.W<-mean(rare.weights)
      }
    }
    rare.Mpvals<-order_Burden(rare.G,obj,xc, weights.beta = Rweights.beta, weights = rare.weights)
    
  }else{
    rare.W=NULL
    rare.Mpvals=NULL
  }
  
  ### common-variant p-values and weights
  if (length(common.G)==0) {
    Mpvals=rare.Mpvals
    W=rare.W
    
  } else {
    if (is.null(common.weights)){
      if (is.null(ncol(common.G))){
        common.MAF=mean(common.G)
      }else {
        common.MAF=colMeans(common.G)/2
      }
      
      common.W<-(dbeta(common.MAF,Cweights.beta[1],Cweights.beta[2])*sqrt(common.MAF*(1-common.MAF)))^2
      
    }else{
      common.W<-common.weights
      
    }
    common.Mpvals<-ord_data_pv(ord_y = obj,x1=common.G,xc=xc)
    
    Mpvals=c(common.Mpvals,rare.Mpvals)
    W=c(common.W,rare.W)
  }
  
  
  ###ACAT
  is.keep<-rep(T,length(Mpvals))
  is.keep[which(Mpvals==1)]<-F  ## remove p-values of 1.
  pval<-ACAT(Mpvals[is.keep],W[is.keep])
 
  return(pval)
}
###################
ord_data_pv=function(ord_y,x1,xc){
  ord_data <- data.frame(ord_y,x1,xc )
  fit0 <- polr(as.ordered(ord_y) ~ 1+xc[,1]+xc[,2], Hess=T, data=ord_data,method = "probit") #H0 
  if (is.null(ncol(x1))){
    fit1 <- polr(as.ordered(ord_y) ~ x1+xc[,1]+xc[,2], Hess=T, data=ord_data,method = "probit") #H1
    A=anova(fit1,fit0)
    pval=A$Pr[2] 
  }else {
    pval=matrix(data = NA,nrow = ncol(x1),ncol = 1)
    
    for(i in 1:ncol(x1))
    {
      fit1 <- polr(as.ordered(ord_y) ~ x1[,i]+xc[,1]+xc[,2], Hess=T, data=ord_data,method = "probit") #H1
      A=anova(fit1,fit0)
      pval[i]=A$Pr[2] 
    }
   
  }
  
  pval
  
}
###################
order_Burden<-function(G,obj,xc,kernel="linear.weighted",weights.beta=c(1,25),weights=NULL){
  
  ### MAF
  MAF<-Matrix::colSums(G)/(2*dim(G)[1])
  p<-length(MAF)
  #### weights
  if (kernel=="linear.weighted"){
    if (is.null(weights)){
      W<-dbeta(MAF,weights.beta[1],weights.beta[2])
    }else{
      if (length(weights)==p){
        W<-weights
      }else{
        stop("The length of weights must equal to the number of variants!")
      }
    }
    
  }else if (kernel=="linear"){
    W<-rep(1,p)  #equal weight
  }else{
    stop("The kernel name is not valid!")
  }
  
  ###### if G is sparse or not
  if (class(G)=="matrix" || class(G)=="dgCMatrix"){
    Gw<-G%*%W
  }else{
    stop("The class of G must be matrix or dgCMatrix!")
  }
  
  ord_data <- data.frame(G,xc, obj)
  fit0 <- polr(as.ordered(obj) ~ 1+xc[,1]+xc[,2], Hess=T, data=ord_data,method = "probit") #H0 
  fit1 <- polr(as.ordered(obj) ~ Gw+xc[,1]+xc[,2], Hess=T, data=ord_data,method = "probit") #H1
  #H1 
  A=anova(fit1,fit0) 
  #B=str(A)
  pval=A$Pr[2] 
  return(pval)
}

