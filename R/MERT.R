## This function conducts the maximin efficiency robust test (MERT) to a given case-control table. ##
## The MERT achieves the maximum minimum efficiency. ##
## Under some conditions, the MERT can be written as the weighted average of two normally distributed tests with the minimum correlation. ##
## In case-control association studies, the extreme pair corresponds to the CATTs under the REC and DOM models. ##
## Under the null hypothesis of no association, the MERT follows the standard normal distribution N(0,1). ##
## The MERT calculates the statistic and associated p-value as well as reporting the conclusion of the hypothesis test. ## 
## data is a 2 by 3 contingency table for analysis. ##
## the rows of the data represent disease status and the columns represent genotypes. ##

MERT <-
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
  ## Define the CATT ##
  CATT=function(data,x=0.5){
    score=c(0,x,1)    
    ## Calculate the test statistic. ## 
    nr=apply(data,2,sum)
    n=sum(nr)
    Rbar=sum(nr*score)/n
    s2=sum(nr*(score-Rbar)^2)
    phi=sum(data[1,])/n
    catt_v=sum(data[1,]*(score-Rbar))/sqrt(phi*(1-phi)*s2)
    return(catt_v)
  }
  ## Calculate the correlation ##
  n=apply(data,2,sum)
  nn=sum(n)
  rho=sqrt((n[1]*n[3])/((nn-n[1])*(nn-n[3])))
  ## Caculate the test statistic of MERT. ##
  catt0=CATT(data,0)
  catt1=CATT(data,1)
  re1=(catt0+catt1)/sqrt(2*(1+rho))
  ## Caculate the p-value. ##
  re2=1-pchisq(re1^2,df=1)
  ## Print the output. ##
  names(re1)="statistic"
  structure(list(statistic=re1,p.value=re2,method="The maximin efficiency robust test",data.name=DNAME),class="htest")

}

## Example 
## ca=c(139,249,112) 
## co=c(136,244,120) 
## a=rbind(ca,co) 
## MERT(a) 
## The maximin efficiency robust test
## data:  a 
## statistic = -0.4962, p-value = 0.6198
