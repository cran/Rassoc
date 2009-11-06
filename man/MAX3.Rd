\name{MAX3}
\Rdversion{1.1}
\alias{MAX3}
\title{
The MAX3 test
}
\description{
This function conducts the MAX3 test to a 2 by 3 case-control contingency table using either empirical or asymptotic methods. This test takes the maximum of the absolute values of the Cochran-Armitage trend tests optimal for recessive, additive/multiplicative and dominant models respectively. 
}
\usage{
MAX3(data, method, m)
}
\arguments{
  \item{data}{
data is a 2 by 3 case-control contingency table. The first and second rows represent the case group and control group respectively. The first, second and third columns represent the genotypes of a susceptibility diallelic marker containing 0, 1 and 2 risk allele respectively. Thus, the numbers in this table represent the genotype counts belonging to the corresponding genotypes and case-control status. 
}
  \item{method}{
method is the approach used for calculating the p-value of MAX3 statistics. boot and bvn are empirical methods and asy is asymptotic method.  
}
  \item{m}{
m is the replication times used for calculating the empirical p-values (boot and bvn). It can be any positive integer (1, for example) while choosing asymptotic method (asy). 
}
}
\details{
boot method re-samples case-control data under null hypothesis and applies MAX3 to each replicate. bvn directly generates statistics of MAX3 rather than case-control data under null hypothesis. asy calculates the p-value from the asymptotic distribution of MAX3 under null hypothesis. 
}
\value{
\item{method}{

the method chosen to calculate the p-value of MAX3}
\item{statistic}{

the statistic of the MAX3 test}
\item{p.value}{

the associated p-value of the MAX3 test}
}
\references{
Freidlin, B, Zheng, G, Li, Z and Gastwirth, JL (2002). Trend tests for case-control studies of genetic markers: power, sample size and robustness. Human Heredity 53, 146-152.

Zang Y, Fung WK and Zheng G(2009). Simple Algorithms to calculate the asymptotic null distribution for MAX3 and genetic model selection in case-control genetic association studies in R. Journal of Statistical software (revised).
}
\examples{
library(Rassoc)
ca=c(139,249,112) 
co=c(136,244,120) 
a=rbind(ca,co) 

MAX3(a,"boot",100000)
## The MAX3 test using the boot method
## data:  a 
## statistic = 0.5993, p-value = 0.7936

MAX3(a,"bvn",100000)
## The MAX3 test using the bvn method
## data:  a 
## statistic = 0.5993, p-value = 0.792

MAX3(a,"asy",1)
## The MAX3 test using the asy method
## data:  a 
## statistic = 0.5993, p-value = 0.7933
}
\seealso{\code{\link{CATT}}}
         

