ABT <-
function(data){
  if((data[1,1]<0.5)||(data[1,2]<0.5)||(data[1,3]<0.5)||(data[2,1]<0.5)||(data[2,2]<0.5)||(data[2,3]<0.5)){
    data=data+matrix(rep(0.5,6),nrow=2)
  }
  rr=data[1,]
  ss=data[2,]
  nn=apply(data,2,sum)
  r=sum(rr)
  s=sum(ss)
  n=sum(nn)
  pd=(2*rr[3]+rr[2])/(2*r)
  ph=(2*ss[3]+ss[2])/(2*s)
  p=(2*nn[3]+nn[2])/(2*n)
  u=pd-ph
  v=p*(1-p)*(1/(2*r)+1/(2*s))
  re1=u/sqrt(v)
  re2=1-pchisq(re1^2,df=1)
  if(re2>=0.05){
    re3="null hypothesis: association doesn't exist under significant level 0.05"
  }
  if(re2<0.05){
    re3="alternative hypothesis: association exists under significant level 0.05"
  }
  out=list("statistics"=re1, "Pvalue"=re2, "conclusion"=re3)
  print("the allelic based test")
  return(out)
}

