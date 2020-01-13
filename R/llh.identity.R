# Log liklihood function: identity link  

llh.identity=function(par,x,w){
  nonzero=ifelse(x!=0,1,x)
  K=length(par)
  
  if(sum(nonzero)==length(nonzero)){
    value=sum(x/(w%*%par[-K])+log(w%*%par[-K]))
  }else{
    
    value=sum((nonzero-1)*log(1-par[K])+par[K]*x/(w%*%par[-K])+nonzero*(log(w%*%par[-K])-2*log(par[K])))
  }
    
  return(value)  
}