#' Group and covariate effects on lncRNA counts by Generalized Linear Model
#'
#' ZIQML.fit estimates the group effect on gene expression using zero-inflated exponential quasi likelihood.  
#' @param edata Normalized counts matrix with genes in rows and samples in columns.
#' @param design.matrix Design matrix for groups and covariates, generated by model.matrix().
#' @param link Link function for the generalized linear model and likelihood function,either 'log' or 'identity'. The default is 'log'. 
#' @return 
#'    \item{Estimates}{Estimated group effect on gene expression by zero-inflated exponential quasi maximum likelihood (ZIQML) estimator.}
#'    \item{logLikelihood}{The value of zero-inflated quasi likelihood.} 
#'    \item{edata}{lncRNA counts or expression matrix.} 
#'    \item{design.matrix}{The design matrix of groups and covariates.} 
#'    \item{link}{The specified link function.} 
#' @examples 
#' 
#' data('hnsc.edata','design') 
#' # 'hnsc.edata' contains FPKM of 1000 lncRNA genes and 80 samples 
#' # 'design' is the design matrix for tissue and batch.
#' 
#' # For the first 100 genes  
#' # Fit GLM by ZIQML with logarithmic link function                                      
#' fit.log=ZIQML.fit(edata=hnsc.edata[1:100,],design.matrix=design,link='log') 
#' 
#' # Fit GLM by ZIQML with identity link function
#' fit.identity=ZIQML.fit(edata=hnsc.edata[1:100,],design.matrix=design,link='identity') 
#' 
#' @export 
#' @importFrom stats optim p.adjust pchisq aggregate as.formula model.matrix

ZIQML.fit=function(edata,design.matrix,link='log'){
  
  if(any(is.na(edata))) stop('Invalid values in expression data: NA')
  if(any(edata<0)) stop('Invalid values in expression data: Negative')
  if(any(is.na(design.matrix))) stop('Design matrix contains NA')
  if(!is.matrix(design.matrix)) stop('Design matrix not in correct form')
  if(nrow(design.matrix)!=ncol(edata)) stop('Dimensions of input matrices do not match')
  if(sum(design.matrix[,1]==0)>0)stop('Design matrix must contain intercept')
  if(!link %in% c('log','identity'))stop('Link function is not valid')
  
  est=exact.llh.value=NULL
  edata=as.matrix(edata)
  design.matrix=as.matrix(design.matrix)
  for(i in 1:nrow(edata)){
    ini=c(mean(edata[i,]),rep(0,(ncol(design.matrix)-1)))
    
    if(link=='identity'){
    mfit=optim(par = ini,
               fn = function(par)llh.identity(c(par,sum(edata[i,]!=0)/length(edata[i,])),edata[i,],design.matrix))
    }else{ 
      mfit=optim(par = c(log(ini[1]),ini[-1]),
                 fn = function(par)llh.log(c(par,sum(edata[i,]!=0)/length(edata[i,])),edata[i,],design.matrix))
    }
  est=rbind(est,mfit$par)
  exact.llh.value=c(exact.llh.value,-mfit$value)
  }
  rownames(est)=rownames(edata)
  colnames(est)=colnames(design.matrix)
  
  
  output=list(Estimates=est,logLikelihood=exact.llh.value,edata=edata,design.matrix=design.matrix,link=link)
  return(output)  
}

  
