MAX3 <-
function(data,method,m){
  library(mvtnorm)
  p=apply(data,2,sum)/sum(data)
  p0=p[1]
  p1=p[2]
  p2=p[3]
  rho005=p2*(p1+2*p0)/(sqrt(p2*(1-p2))*sqrt((p1+2*p2)*p0+(p1+2*p0)*p2))
  rho105=p0*(p1+2*p2)/(sqrt(p0*(1-p0))*sqrt((p1+2*p2)*p0+(p1+2*p0)*p2))
  rho01=sqrt((p0*p2)/((1-p0)*(1-p2)))
  w0=(rho005-rho01*rho105)/(1-rho01^2)
  w1=(rho105-rho01*rho005)/(1-rho01^2)
  r=sum(data[1,])
  s=sum(data[2,])
  CATT=function(tab,score){
    if((tab[1,1]<0.5)||(tab[1,2]<0.5)||(tab[1,3]<0.5)||(tab[2,1]<0.5)||(tab[2,2]<0.5)||(tab[2,3]<0.5)){
      tab=tab+matrix(rep(0.5,6),nrow=2)
    }
    nr=apply(tab,2,sum)
    n=sum(nr)
    Rbar=sum(nr*score)/n
    s2=sum(nr*(score-Rbar)^2)
    phi=sum(tab[1,])/n
    catt_v=sum(tab[1,]*(score-Rbar))/sqrt(phi*(1-phi)*s2)


    return(catt_v)
  }
  maxr=max(abs(CATT(data,c(0,0,1))),abs(CATT(data,c(0,0.5,1))),abs(CATT(data,c(0,1,1))))
  print("MAX3 test")
  if(method=="boot"){
    CATTN=function(data,score){
      tab=matrix(data,nrow=2,byrow=TRUE)
      if((tab[1,1]<0.5)||(tab[1,2]<0.5)||(tab[1,3]<0.5)||(tab[2,1]<0.5)||(tab[2,2]<0.5)||(tab[2,3]<0.5)){
        tab=tab+matrix(rep(0.5,6),nrow=2)
      }
      nr=apply(tab,2,sum)
      n=sum(nr)
      Rbar=sum(nr*score)/n
      s2=sum(nr*(score-Rbar)^2)
      phi=sum(tab[1,])/n
      catt_v=sum(tab[1,]*(score-Rbar))/sqrt(phi*(1-phi)*s2)


      return(catt_v)
    }
    ca=rmultinom(m,r,p)
    co=rmultinom(m,s,p)
    caco=rbind(ca,co)
    caco0=apply(caco,2,CATTN,score=c(0,0,1))
    caco05=apply(caco,2,CATTN,score=c(0,0.5,1))
    caco1=apply(caco,2,CATTN,score=c(0,1,1))
    cacos=rbind(abs(caco0),abs(caco05),abs(caco1))
    max1=apply(cacos,2,max)
    pv1=length(max1[max1>=maxr])/m
    if(pv1>=0.05){
      conclu="null hypothesis: association doesn't exist under significant level 0.05"
    }
    if(pv1<0.05){
      conclu="alternative hypothesis: association exists under significant level 0.05"
    }
    return(list("method"="boot","statistics"=maxr,"Pvalue"=pv1,"conclusion"=conclu))
  }
  if(method=="bvn"){
    BV05=function(x,a1,a2){
      return(a1*x[1]+a2*x[2])
    }
    bv01=rmvnorm(m,mean=c(0,0),sigma=matrix(c(1,rho01,rho01,1),nrow=2))
    bv05=apply(bv01,1,BV05,a1=w0,a2=w1)
    bv=cbind(abs(bv01),abs(bv05))
    max2=apply(bv,1,max)
    pv2=length(max2[max2>=maxr])/m
    if(pv2>=0.05){
      conclu="null hypothesis: association doesn't exist under significant level 0.05"
    }
    if(pv2<0.05){
      conclu="alternative hypothesis: association exists under significant level 0.05"
    }
    return(list("method"="bvn","statistics"=maxr,"Pvalue"=pv2,"conclusion"=conclu))
  }
  if(method=="asy"){
    fun1=function(t,z){
      return(pnorm((t-rho01*z)/sqrt(1-rho01^2))*dnorm(z))
    }
    fun2=function(t,z){
      return(pnorm(((t-w0*z)/w1-rho01*z)/sqrt(1-rho01^2))*dnorm(z))
    }
    fun3=function(t,z){
      return(pnorm((-t-rho01*z)/sqrt(1-rho01^2))*dnorm(z))
    }
    asy=function(t){
      z1=2*integrate(function(z){fun1(t,z)},lower=0,upper=t*(1-w1)/w0,subdivisions=1000)$value
      z2=2*integrate(function(z){fun2(t,z)},lower=t*(1-w1)/w0,upper=t,subdivisions=1000)$value
      z3=-2*integrate(function(z){fun3(t,z)},lower=0,upper=t,subdivisions=500)$value
      return(z1+z2+z3)
    }
    pv3=1-asy(maxr)
    if(pv3>=0.05){
      conclu="null hypothesis: association doesn't exist under significant level 0.05"
    }
    if(pv3<0.05){
      conclu="alternative hypothesis: association exists under significant level 0.05"
    }
    return(list("method"="asy","statistics"=maxr,"Pvalue"=pv3,"conclusion"=conclu))
  }
}

