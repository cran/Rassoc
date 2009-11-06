## This function conducts the genetic model selection (GMS) to a given case-control table. ##
## The GMS calculates the statistic and associated p-value as well as reporting the conclusion of the hypothesis test. ##
## GMS is a two-phase adaptive test. ##
## In the first phase, the underlying genetic model is selected using the Hardy-Weinberg disequilibrium trend test (HWDTT). ##
## In the second phase, the CATT with the selected genetic model is applied to test for association. ##
## CATT(0.5) is used to detect the risk allele. ##
## data is a 2 by 3 contingency table for analysis. ##
## The rows of the data represent disease status and the columns represent genotypes. ## 
## method is the method for calculating the p-values of MAX3. ##
## "boot" represents the parametric bootstrap method. ##
## "bvn" represents the simulation from the bivariate normal distribution. ##
## "asy" represents the asymptotic null distribution method. ##
## m is the number of replicates for "boot" and "bvn". ##
## m can be any positive integer for "asy". ##
## m MUST be positive intefer ##

GMS <-
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
  ## Define CATT. ##
  CATT=function(tab,score){
    n=sum(tab)
    nr=apply(tab,2,sum)
    Rbar=sum(nr*score)/n
    s2=sum(nr*(score-Rbar)^2)
    phi=sum(tab[1,])/n
    catt_v=sum(tab[1,]*(score-Rbar))/sqrt(phi*(1-phi)*s2)
    ## Report the statistic of CATT. ##
    return(catt_v)
  }
  ## Build the index function. ##
  index=function(a){
    if(a=="TRUE"){
      b=1
    }
    else{
      b=0
    }
    return(b)
  }
  ## Define the HWDTT. ##
  HWDTT=function(tab){
    n=sum(tab)
    r=sum(tab[1,])
    s=sum(tab[2,])
    nn=apply(tab,2,sum)
    pr=tab[1,]/r
    ps=tab[2,]/s
    pn=nn/n
    deltap=pr[3]-(pr[3]+0.5*pr[2])^2
    deltaq=ps[3]-(ps[3]+0.5*ps[2])^2
    u=sqrt(r*s/n)*(deltap-deltaq)
    v1=1-pn[3]-0.5*pn[2]
    v2=pn[3]+0.5*pn[2]
    ## Report the statistic of HWDTT. ##
    return(u/(v1*v2))
  }
  ## Build the GMS function. ##
  REAL=function(data){
    catt0=CATT(data,c(0,0,1))
    catt05=CATT(data,c(0,0.5,1))
    catt1=CATT(data,c(0,1,1))
    se=HWDTT(data)
    c0=qnorm(0.95)
    a1=catt0*index((catt05>0)&&(se>c0))+catt1*index((catt05>0)&&(se<(-c0)))+catt05*index((catt05>0)&&(abs(se)<=c0))
    a2=-catt1*index((catt05<=0)&&(se>c0))-catt0*index((catt05<=0)&&(se<(-c0)))-catt05*index((catt05<=0)&&(abs(se)<=c0))
    return(a1+a2)  
  }
  GMS1=function(data){
    c0=qnorm(0.95)
    catt0=data[1]
    catt1=data[2]
    catt05=data[3]
    se=data[4]
    a1=catt0*index((catt05>0)&&(se>c0))+catt1*index((catt05>0)&&(se<(-c0)))+catt05*index((catt05>0)&&(abs(se)<=c0))
    a2=-catt1*index((catt05<=0)&&(se>c0))-catt0*index((catt05<=0)&&(se<(-c0)))-catt05*index((catt05<=0)&&(abs(se)<=c0))
    return(a1+a2)
  }
  GMS2=function(data){
    c0=qnorm(0.95)
    tab=matrix(data,nrow=2,byrow=TRUE)
    catt0=CATT(tab,c(0,0,1))
    catt05=CATT(tab,c(0,0.5,1))
    catt1=CATT(tab,c(0,1,1))
    se=HWDTT(tab)
    a1=catt0*index((catt05>0)&&(se>c0))+catt1*index((catt05>0)&&(se<(-c0)))+catt05*index((catt05>0)&&(abs(se)<=c0))
    a2=-catt1*index((catt05<=0)&&(se>c0))-catt0*index((catt05<=0)&&(se<(-c0)))-catt05*index((catt05<=0)&&(abs(se)<=c0))
    return(a1+a2)
  }
  ## Caculate the correlations. ##
  MAF=sum(data[,3])/sum(data)+0.5*(sum(data[,2])/sum(data))
  rho01=sqrt(MAF*(1-MAF)/((1+MAF)*(2-MAF)))
  rho005=sqrt(2*(MAF)/(1+MAF))
  rho051=sqrt(2*(1-MAF)/(2-MAF))
  rho0=sqrt((1-MAF)/(1+MAF))
  rho1=-sqrt(MAF/(2-MAF))
  omega0=sqrt(MAF*(1+MAF)/2)
  omega1=sqrt((1-MAF)*(2-MAF)/2)
  phi0=(rho0-rho1*rho01)/(1-rho01^2)
  phi1=(rho1-rho0*rho01)/(1-rho01^2)
  p=c(sum(data[,1]),sum(data[,2]),sum(data[,3]))/sum(data)
  r=sum(data[1,])
  s=sum(data[2,])
  ## Report the statistic of GMS. ##
  real=REAL(data)
  ## Define the asymptotic formula. ##
  ASY=function(t){
    c0=qnorm(0.95)
    s1=matrix(c(1,0,rho005,0,1,rho0,rho005,rho0,1),nrow=3)
    s2=matrix(c(1,0,rho051,0,1,-rho1,rho051,-rho1,1),nrow=3)
    l1=pmvnorm(lower=-Inf,upper=c(0,-c0,-t),mean=c(0,0,0),sigma=s1)[1] 
    l2=pmvnorm(lower=-Inf,upper=c(0,-c0,-t),mean=c(0,0,0),sigma=s2)[1]
    a=2*(l1+l2)+1.8*pnorm(min(0,-t))
    return(a)  
  }
  if((method!="boot")&&(method!="bvn")&&(method!="asy"))
  stop("method must be boot, bvn or asy")
  ## Use "asy" to caculate the p-value of GMS. ##
  if(method=="asy"){
    ## Report asymptotic p-value using "asy" method. ##
    b=ASY(real)
    md="The GMS test using the asy method"
  }
  ## Use "bvn" to caculate the p-value of GMS. ##
  if(method=="bvn"){
    BV05=function(x,a1,a2){
      return(a1*x[1]+a2*x[2])
    }
    bv01=rmvnorm(m,mean=c(0,0),sigma=matrix(c(1,rho01,rho01,1),nrow=2))
    bv05=apply(bv01,1,BV05,a1=omega0,a2=omega1)
    bvh=apply(bv01,1,BV05,a1=phi0,a2=phi1)
    bv=cbind(bv01,bv05,bvh)
    gmstt1=apply(bv,1,GMS1)
    ## Report empirical p-value using "bvn" method. ##
    b=length(gmstt1[gmstt1>real])/m
    md="The GMS test using the bvn method"   
  }
  ## Use "boot" to caculate the p-value of GMS. ##
  if(method=="boot"){
    ca=rmultinom(m,r,p)
    co=rmultinom(m,s,p)
    caco=rbind(ca,co)
    gmstt2=apply(caco,2,GMS2)
    ## Report empirical p-value using "boot" method. ##
    b=length(gmstt2[gmstt2>real])/m
    md="The GMS test using the boot method"    
  }
  ## print the output. ##
  names(real)="statistic"
  structure(list(statistic=real,p.value=b,method=md,data.name=DNAME),class="htest")

}

## Examples 
## ca=c(139,249,112) 
## co=c(136,244,120) 
## a=rbind(ca,co) 
## GMS(a,"boot",100000.5) 
## Error in GMS(a, "boot", 100000.5) : warning: m MUST be integer 
## GMS(a,"hansi",10000)
## Error in GMS(a, "hansi", 10000) : method must be boot, bvn or asy

## GMS(a,"boot",100000) 
## The GMS test using the boot method
## data:  a 
## statistic = 0.4894, p-value = 0.6658

## GMS(a,"bvn",100000)
## The GMS test using the bvn method
## data:  a 
## statistic = 0.4894, p-value = 0.6598

## GMS(a,"asy",100000)
## The GMS test using the asy method
## data:  a 
## statistic = 0.4894, p-value = 0.6621



