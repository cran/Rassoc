## This function conducts the MAX3 to a given case-control table. ##
## MAX3 takes the maximum of the absolute values of the three CATTs respectively optimal for the REC, ADD and DOM models. ##
## The MAX3 calculates the statistic and associated p-value as well as reporting the conclusion of the hypothesis test. ##
## data is a 2 by 3 contingency table for analysis. ##
## The rows of the data represent disease status and the columns represent genotypes. ## 
## method is the method for calculating the p-values of MAX3. ##
## "boot" represents the parametric bootstrap method. ##
## "bvn" represents the simulation from the bivariate normal distribution. ##
## "asy" represents the asymptotic null distribution method. ##
## m is the number of replicates for "boot" and "bvn". ##
## m can be any positive integer for "asy". ##
## m MUST be positive intefer ##

MAX3 <-
function(data,method="asy",m=1){
  DNAME=deparse(substitute(data))
  data=as.matrix(data)
  ## Check if data is a 2 by 3 contingency table. ##
  datanew=table(data,exclude=c(NA, NaN, Inf, -Inf))
  if(sum(datanew)!=6)
  stop("data must be a 2 by 3 table without infinite and missing values.")
  if(any(data<0))
  stop("all entries of data must be non-negative.")
  if((dim(data)[1]!=2)||(dim(data)[2]!=3))
  stop("data must be a 2 by 3 table.")
  if((data[1,1]<0.5)||(data[1,2]<0.5)||(data[1,3]<0.5)||(data[2,1]<0.5)||(data[2,2]<0.5)||(data[2,3]<0.5)){
  warning("At least one cell of the table is zero.")
  data=data+matrix(rep(0.5,6),nrow=2)
  }
  ## Check if m is integer. ##
  if ((trunc(m)-m)<0)
  stop("m MUST be integer")
  ## Caculate the correlations. ##
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
  ## Define CATT. ##
  CATT=function(tab,score){
    nr=apply(tab,2,sum)
    n=sum(nr)
    Rbar=sum(nr*score)/n
    s2=sum(nr*(score-Rbar)^2)
    phi=sum(tab[1,])/n
    catt_v=sum(tab[1,]*(score-Rbar))/sqrt(phi*(1-phi)*s2)
    ## Report the statistic of CATT. ##
    return(catt_v)
  }
  ## Report the statistic of MAX3. ##
  maxr=max(abs(CATT(data,c(0,0,1))),abs(CATT(data,c(0,0.5,1))),abs(CATT(data,c(0,1,1))))
  if((method!="boot")&&(method!="bvn")&&(method!="asy"))
  stop("method must be boot, bvn or asy.")
  ## Use "boot" to caculate the p-value of MAX3. ##
  if(method=="boot"){
    CATTN=function(data,score){
      tab=matrix(data,nrow=2,byrow=TRUE)
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
    ## Report empirical p-value using "boot" method. ##
    pv=length(max1[max1>=maxr])/m
    md="The MAX3 test using the boot method"
   
  }
  ## Use "bvn" to caculate the p-value of MAX3. ##
  if(method=="bvn"){
    BV05=function(x,a1,a2){
      return(a1*x[1]+a2*x[2])
    }
    bv01=rmvnorm(m,mean=c(0,0),sigma=matrix(c(1,rho01,rho01,1),nrow=2))
    bv05=apply(bv01,1,BV05,a1=w0,a2=w1)
    bv=cbind(abs(bv01),abs(bv05))
    max2=apply(bv,1,max)
    ## Report empirical p-value using "bvn" method. ##
    pv=length(max2[max2>=maxr])/m
    md="The MAX3 test using the bvn method"
    
  }
  ## Use "asy" to caculate the p-value of MAX3. ##
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
    ## Report asymptotic p-value using "asy" method. ##
    pv=1-asy(maxr)
    md="The MAX3 test using the asy method"
    
  }
  ## print the output. ##
    names(maxr)="statistic"
    structure(list(statistic=maxr,p.value=pv,method=md,data.name=DNAME),class="htest")

}

## Examples 
## ca=c(139,249,112) 
## co=c(136,244,120) 
## a=rbind(ca,co) 
## MAX3(a,"boot",100000.5) 
## Error in MAX3(a, "boot", 100000.5) : warning: m MUST be integer 
## MAX3(a,"hansi",100000) 
## Error in MAX3(a, "hansi", 1e+05) : method must be boot, bvn or asy.
## MAX3(a,"boot",100000)
## The MAX3 test using the boot method
## data:  a 
## statistic = 0.5993, p-value = 0.7936
## MAX3(a,"bvn",100000)
## The MAX3 test using the bvn method
## data:  a 
## statistic = 0.5993, p-value = 0.792
## MAX3(a,"asy",1)
## The MAX3 test using the asy method
## data:  a 
## statistic = 0.5993, p-value = 0.7933
