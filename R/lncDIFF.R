#' lncRNA Differential Expression (DE) analysis  
#'
#' lncDIFF returns DE analysis results based on lncRNA counts and grouping variables.  
#' @param edata Normalized counts matrix with genes in rows and samples in columns.
#' @param group Primary factor of interest in DE analysis, e.g., treatment groups, tissue types, other phenotypes. 
#' @param covariate Other variables (or covariates) associated with expression level. Input must be a matrix or data frame with each column being a covariate matching to \code{group}  
#' @param link.function Link function for the generalized linear model, either 'log' or 'identity', default as 'log'.
#' @param CompareGroups Labels of treatment groups or phenotypes of interest to be compared in DE analysis. Input must be a vector of \code{group} labels without duplicates. 
#' @param simulated.pvalue If empirical p-values are computed, simulated.pvalue=TRUE. The default is FALSE. 
#' @param permutation The number of permutations used in simulating pvalues. The default value is 100. 
#' @return 
#'        \item{DE.results}{Likelihood ratio test results with test statistics, p-value, FDR, DE genes, groupwise mean expression, fold change (if two groups are compared). 
#'        If simulated.pvalue=TRUE, test.results also includes simulated p-value and FDR.}  
#'        \item{full.model.fit}{Generalized linear model with zero-inflated Exponential likelihood function, estimating group effect compared to a reference group.}
#' @references Li, Q., Yu, X., Chaudhary, R. et al.'lncDIFF: a novel quasi-likelihood method for differential expression analysis of non-coding RNA'. BMC Genomics (2019) 20: 539.     
#' @examples
#' 
#' data('hnsc.edata','tissue','cov')  
#' 
#' # DE analysis comparing two groups (normal vs tumor) for 100 genes
#' result=lncDIFF(edata=hnsc.edata[1:100,],group=tissue,covariate=cov) 
#' 
#' # Recommend at least 50 permutations if simulated.pvalue=TRUE
#' 
#' 
#' 
#' 
#' @export  
#' @importFrom stats optim p.adjust pchisq aggregate as.formula model.matrix

lncDIFF<-function(edata,group,covariate=NULL,link.function='log',CompareGroups=NULL,simulated.pvalue=FALSE,permutation=100){
  group.labels<-names(table(group))
  if(is.null(CompareGroups)) {
    cat('Compared groups are not specified, default to all groups','\n')
    CompareGroups<-group.labels
  }  
  if(length(CompareGroups)>length(group.labels)) stop('Duplicates or unspecified groups in CompareGroups')
  if(sum(!CompareGroups %in% group)>0) stop('Groups to be compared (CompareGroups) are not in the range of labels')
  if(length(CompareGroups)==1) stop('Specify at least 2 groups to be compared')
  if(!is.vector(group)) stop('Treatment groups or phenotypes of interest (group) must be a vector ')
  if(!is.null(covariate)){ 
    if(length(group)!=nrow(covariate)) 
    stop('Dimensions of covariate and group do not match')
  }
  if(length(group)!=ncol(edata)) stop('Dimensions of counts (edata) and group do not match')
  
  
  
  if(length(CompareGroups)>2) cat('More than 2 groups are compared, fold change are not computed','\n')
  
  full.coefficients=which(group.labels %in% CompareGroups)
  
  if(! 1 %in% full.coefficients){
    group<-factor(group,levels = c(CompareGroups,group.labels[-full.coefficients]))
    group.labels<-names(table(group))
    full.coefficients=which(group.labels %in% CompareGroups)
  }
  
  pdata=as.data.frame(cbind(group,covariate))
  formula='~'
  for(i in colnames(pdata)[-ncol(pdata)]){
      formula=paste(formula,i,'+',sep = '')
  }
  formula=paste(formula,colnames(pdata)[ncol(pdata)],sep='')
  
  design.matrix=model.matrix(as.formula(formula),pdata)
  colnames(design.matrix)[2:length(group.labels)]=paste(deparse(substitute(group)),group.labels[-1],sep = '')
  ZIQML.fit.full=ZIQML.fit(edata,design.matrix,link=link.function)
  
  test.coef=sort(full.coefficients)[-1]
  
  n=nrow(design.matrix)
  g=nrow(edata)
  test=LRT(ZIQML.fit.full,coef=test.coef)
  LRT.stat=test$LRT.stat
  LRT.pvalue=test$LRT.pvalue
  LRT.fdr=p.adjust(LRT.pvalue,method = 'BH')
  DE.Gene=ifelse(LRT.fdr<0.05,'Yes','No')
  results=data.frame(test.statistics=LRT.stat,Pvalue=LRT.pvalue,FDR=LRT.fdr,DE.Gene=DE.Gene)
  
  compare.id=group %in% CompareGroups
  sub.edata=t(edata[,compare.id])
  groupwise.mean=aggregate(sub.edata,by=list(group[compare.id]),FUN=mean)
  rownames(groupwise.mean)=groupwise.mean[,1]
  
  
  groupwise.mean=t(groupwise.mean[,-1])
  colnames(groupwise.mean)=paste('Mean',colnames(groupwise.mean),sep ='_' )
  results=cbind(results,groupwise.mean)
  
  if(simulated.pvalue){
    LRT.STAT=NULL
    for(i in 1:permutation){
      id=sample(1:n,n)
      ZIQML.fit.null=ZIQML.fit(edata,design.matrix[id,],link = ZIQML.fit.full$link)
      test=LRT(ZIQML.fit.null,coef=test.coef)
      LRT.STAT=cbind(LRT.STAT,test$LRT.stat)
      
    }
    
    LRT.simulated.pvalue=(0.1+rowSums(LRT.STAT>LRT.stat))/permutation
    LRT.simulated.pvalue=lapply(LRT.simulated.pvalue,function(x)min(x,1))
    LRT.simulated.fdr=p.adjust(LRT.simulated.pvalue,method = 'BH')
    results$Simulated.Pvalue=LRT.simulated.pvalue
    results$Simulated.FDR=LRT.simulated.fdr
    results$DE.Gene.Simulated.Fdr=ifelse(results$Simulated.FDR<0.05,'Yes','No')
    
  } 
  
  
  if(length(CompareGroups)==2){
    results$Fold.Change=groupwise.mean[,1]/groupwise.mean[,2]
    results$Log2.Fold.Change=log2(groupwise.mean[,1]/groupwise.mean[,2])  
  }
  
  rownames(results)=rownames(ZIQML.fit.full$edata)
  output=list(DE.results=results,full.model.fit=ZIQML.fit.full)
  return(output)
}  