## This function conducts the Cochran-Armitage trend test (CATT) to a given case-control table. ## 
## The CATT detects association by comparing the genotype frequencies between cases and controls. ##
## Under the null hypothesis of no association, the CATT follows the standard normal distribution N(0,1). ##
## The CATT calculates the statistic and associated p-value as well as reporting the conclusion of the hypothesis test. ## 
## data is a 2 by 3 contingency table for analysis. ##
## the rows of the data represent disease status and the columns represent genotypes. ## 
## x is the score for the CATT. It can be any real number between 0 and 1. ##
## Specifically, x=0, 0.5 and 1 are optimal for REC, MUL/ADD and DOM models respectively. ##

CATT <-
function(data,x=0.5){
  score=c(0,x,1)
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
  ## Calculate the test statistic. ## 
  nr=apply(data,2,sum)
  n=sum(nr)
  Rbar=sum(nr*score)/n
  s2=sum(nr*(score-Rbar)^2)
  phi=sum(data[1,])/n
  re1=sum(data[1,]*(score-Rbar))/sqrt(phi*(1-phi)*s2)
  ## Calculate the p-value. ##
  re2=1-pchisq(re1^2,df=1)
  ## Print the output. ##
  names(re1)="statistic"
  structure(list(statistic=re1,p.value=re2,method="The Cochran-Armitage trend test",data.name=DNAME),class="htest")

}

## Example 
## ca=c(139,249,112) 
## co=c(136,244,120) 
## a=rbind(ca,co) 
## CATT(a,0.5) 
## the Cochran-Armitage trend test
## data:  a 
## statistic = -0.4894, p-value = 0.6245
