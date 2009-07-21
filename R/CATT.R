CATT <-
function(data,x){
  score=c(0,x,1)
  if((data[1,1]<0.5)||(data[1,2]<0.5)||(data[1,3]<0.5)||(data[2,1]<0.5)||(data[2,2]<0.5)||(data[2,3]<0.5)){
    data=data+matrix(rep(0.5,6),nrow=2)
  }
  nr=apply(data,2,sum)
  n=sum(nr)
  Rbar=sum(nr*score)/n
  s2=sum(nr*(score-Rbar)^2)
  phi=sum(data[1,])/n
  re1=sum(data[1,]*(score-Rbar))/sqrt(phi*(1-phi)*s2)
  re2=1-pchisq(re1^2,df=1)
  if(re2>=0.05){
    re3="null hypothesis: association doesn't exist under significant level 0.05"
  }
  if(re2<0.05){
    re3="alternative hypothesis: association exists under significant level 0.05"
  }
  out=list("statistics"=re1, "Pvalue"=re2, "conclusion"=re3)
  print("the Cochran-Armitage trend test")
  return(out)
}

