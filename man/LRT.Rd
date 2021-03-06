% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LRT.R
\name{LRT}
\alias{LRT}
\title{Likelihood ratio test based on ZIQML.fit()}
\usage{
LRT(ZIQML.fit, coef = NULL)
}
\arguments{
\item{ZIQML.fit}{Object returned by ZIQML.fit()}

\item{coef}{An integer or vector indicating the coefficient(s) in design matrix to be tested. coef=1 is the intercept (i.e. baseline group effect),
and should not be tested.}
}
\value{
\item{LRT.stat}{Likelihood ratio test statistics.}  
       \item{LRT.pvalue}{Likelihood ratio test p-value.}
}
\description{
ZIQML.LRT returns the likelihood ratio test statistics and p-value based on the object returned by ZIQML.fit().
}
\examples{

data('hnsc.edata','design')  
# 'hnsc.edata' contains FPKM of 1132 lncRNA genes and 80 samples. 
# 'design' is the design matrix of tissue type (tumor vs normal). 

# Fit GLM by ZIQML.fit for the first 100 genes 
fit.log=ZIQML.fit(edata=hnsc.edata[1:100,],design.matrix=design) 


# Likelihood ratio test to compare tumor vs normal in gene expression level. 
LRT.results=LRT(fit.log,coef=2)  

}
