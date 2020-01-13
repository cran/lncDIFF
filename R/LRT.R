#' Likelihood ratio test based on ZIQML.fit()
#'
#' ZIQML.LRT returns the likelihood ratio test statistics and p-value based on the object returned by ZIQML.fit().
#' 
#' @param ZIQML.fit Object returned by ZIQML.fit()
#' @param coef An integer or vector indicating the coefficient(s) in design matrix to be tested. coef=1 is the intercept (i.e. baseline group effect),
#' and should not be tested. 
#' 
#' @return 
#'        \item{LRT.stat}{Likelihood ratio test statistics.}  
#'        \item{LRT.pvalue}{Likelihood ratio test p-value.}
#'      
#' @examples
#' 
#' data('hnsc.edata','design')  
#' # 'hnsc.edata' contains FPKM of 1132 lncRNA genes and 80 samples. 
#' # 'design' is the design matrix of tissue type (tumor vs normal). 
#' 
#' # Fit GLM by ZIQML.fit for the first 100 genes 
#' fit.log=ZIQML.fit(edata=hnsc.edata[1:100,],design.matrix=design) 
#' 
#' 
#' # Likelihood ratio test to compare tumor vs normal in gene expression level. 
#' LRT.results=LRT(fit.log,coef=2)  
#' 
#' @export 
#' @importFrom stats optim p.adjust pchisq aggregate as.formula model.matrix

LRT<-function(ZIQML.fit,coef=NULL){
  m=length(coef)
  design.matrix=ZIQML.fit$design.matrix
  edata=ZIQML.fit$edata
  LRT.stat= LRT.pvalue= LRT.fdr=NULL
  for(i in 1:nrow(edata)){
    ini=c(mean(edata[i,]),rep(0,(ncol(design.matrix)-1)))
    ini.alt=ini
    ini.alt=ini.alt[-coef]
    
    if(ZIQML.fit$link=='identity'){
      
      if(length(ini.alt)==1){
        reduce=optim(par = ini.alt,fn = function(y){
          x=1:ncol(design.matrix)
          x[coef]=0
          x[-coef]=y
          llh.identity(c(x,sum(edata[i,]!=0)/length(edata[i,])),edata[i,],design.matrix)},method = 'Brent',
          lower=-10^5,upper=10^5)
        
      }else{
        reduce=optim(par = ini.alt,fn = function(y){
          x=1:ncol(design.matrix)
          x[coef]=0
          x[-coef]=y
          llh.identity(c(x,sum(edata[i,]!=0)/length(edata[i,])),edata[i,],design.matrix)})
      }
    }else{
      if(length(ini.alt)==1){
        
        reduce=optim(par = c(log(ini.alt[1]),ini.alt[-1]),fn = function(y){
          x=1:ncol(design.matrix)
          x[coef]=0
          x[-coef]=y
          llh.log(c(x,sum(edata[i,]!=0)/length(edata[i,])),edata[i,],design.matrix)},method = 'Brent',
          lower=log(10^(-5)),upper=log(10^5))
        
      }else{
        
        reduce=optim(par = c(log(ini.alt[1]),ini.alt[-1]),fn = function(y){
          x=1:ncol(design.matrix)
          x[coef]=0
          x[-coef]=y
          llh.log(c(x,sum(edata[i,]!=0)/length(edata[i,])),edata[i,],design.matrix)})
      }
    }
    
    lr.par=2*(ZIQML.fit$logLikelihood[i]+reduce$value)
    
    lr.p=pchisq(lr.par,1,lower.tail = F)
    
    LRT.stat=c(LRT.stat,lr.par)
    LRT.pvalue=c(LRT.pvalue,lr.p)
  }
    output=list(LRT.stat=LRT.stat,LRT.pvalue=LRT.pvalue)
    return(output)

}