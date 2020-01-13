# Log liklihood function: logarithmic link  

llh.log=function(par,x,w){
  nonzero=ifelse(x!=0,1,x)
  K=length(par)
  if(sum(nonzero)==length(nonzero)){
    value=sum(x*exp(-(w%*%par[-K]))+w%*%par[-K])
  }else{
    
    value=sum((nonzero-1)*log(1-par[K])+par[K]*x*exp(-(w%*%par[-K]))+nonzero*(w%*%par[-K]-2*log(par[K])))
  }
    return(value)  
}