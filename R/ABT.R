## This function conducts the allelic based test (ABT) to a given case-control table. ##
## The ABT detects association by comparing the allele frequencies between cases and controls. ##
## Under the null hypothesis of no association, the ABT follows the standard normal distribution N(0,1). ##
## The ABT needs the Hardy-Weinberg Equilibrium (HWE) proportions to hold in the population. ##  
## The ABT calculates the statistic and associated p-value as well as reporting the conclusion of the hypothesis test. ## 
## data is a 2 by 3 contingency table for analysis. ##
## the rows of the data represent disease status and the columns represent genotypes. ## 

ABT <-
function(data){
  DNAME=deparse(substitute(data))
  ## Check if data is a 2 by 3 contingency table. ##
  if(any(is.na(data)==TRUE)||any(abs(data)>10^9))  
  stop("data must be a 2 by 3 table without infinite and missing values.")
  if(any(data<0))
  stop("all entries of data must be non-negative.")
  if((dim(data)[1]!=2)||(dim(data)[2]!=3))
  stop("data must be a 2 by 3 table.")
  if((data[1,1]<0.5)||(data[1,2]<0.5)||(data[1,3]<0.5)||(data[2,1]<0.5)||(data[2,2]<0.5)||(data[2,3]<0.5)){
    warning("At least one cell of the table is zero.")
    data=data+matrix(rep(0.5,6),nrow=2)
  }
  ## rr is the vector of case group. ##
  rr=data[1,]
  ## ss is the vector of control group. ##
  ss=data[2,]
  ## nn is the vector of samples. ##
  nn=apply(data,2,sum)
  ## r is the number of cases. ##
  r=sum(rr)
  ## s is the number of controls. ##
  s=sum(ss)
  ## n is the number of samples. ##
  n=sum(nn)
  ## pd is the estimated allele frequency in cases. ##
  pd=(2*rr[3]+rr[2])/(2*r)
  ## ph is the estimated allele frequency in controls. ##
  ph=(2*ss[3]+ss[2])/(2*s)
  ## p is the estimated allele frequency under null hypothesis. ##
  p=(2*nn[3]+nn[2])/(2*n)
  ## Calculate the test statistic. ##
  u=pd-ph
  v=p*(1-p)*(1/(2*r)+1/(2*s))
  re1=u/sqrt(v)
  ## Calculate the p-value. ##
  re2=1-pchisq(re1^2,df=1)
  ## Print the output. ##
  names(re1)="statistic"
  structure(list(statistic=re1,p.value=re2,method="The allelic based test",data.name=DNAME),class="htest")
  
}



## Example 
## ca=c(139,249,112) 
## co=c(136,244,120) 
## a=rbind(ca,co) 
## ABT(a) 
## "The allelic based test" 
## data: a 
## statistic=-0.4924, p-value=0.6224



