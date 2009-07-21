\name{MAX3}
\Rdversion{1.1}
\alias{MAX3}
\title{
The MAX3 test
}
\description{
This function conducts the MAX3 test to a 2 by 3 case-control contingency table using either empirical or asymptotic methods. This test takes the maximum of the absolute values of the Cochran-Armitage trend tests optimal for recessive, additive/multiplicative and dominant models respectively. The method chosen, test statistics, associated p-value and conclusion for the hypothesis test are listed as results.   
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
\item{statistics}{

the statistics of the MAX3 test}
\item{Pvalue}{

the associated p-value of the MAX3 test}
\item{conclusion}{

the conclusion of the hypothesis test}
}
\references{
Freidlin, B, Zheng, G, Li, Z and Gastwirth, JL (2002). Trend tests for case-control studies of genetic markers: power, sample size and robustness. Human Heredity 53, 146-152.

Zang Y, Fung WK and Zheng G(2009). Simple Algorithms to calculate the asymptotic null distribution for MAX3 and genetic model selection in case-control genetic association studies in R. Journal of Statistical software (revised).
}
\examples{
library(Rassoc)
data(caco)
ex=matrix(caco[1,],nrow=2,byrow=TRUE)
## ex is an example of a 2 by 3 case-control contingency table. 

## Using boot method 
MAX3(ex,"boot",10000)
## Conduct the MAX3 test to ex, boot method is used based on 1000000   
##replicates to approximate the p-value of MAX3
## "MAX3 test"
##$method
## "boot"
## the method chosen is boot
##$statistics
## 5.117125
## The statistics of the test is 5.117125
##$Pvalue
## 1e-06
## The associated p-value of the test using boot method is 1e-06
##$conclusion
## "alternative hypothesis: association exists under significant level 0.05"
## Under the significant level 0.05, the alternative hypothesis that the 

##association exists is accepted according to the given p-values.

## Using bvn method 
MAX3(ex,"bvn",10000)
## Conduct the MAX3 test to ex, bvn method is used based on 1000000   
##replicates to approximate the p-value of MAX3
## "MAX3 test"
##$method
## "bvn"
## the method chosen is bvn
##$statistics
## 5.117125
## The statistics of the test is 5.117125
##$Pvalue
## 3e-06
## The associated p-value of the test using bvn method is 3e-06
##$conclusion
## "alternative hypothesis: association exists under significant level 0.05"
## Under the significant level 0.05, the alternative hypothesis that the 

##association exists is accepted according to the given p-values.

## Using asy method 
MAX3(ex,"asy",1)
## Conduct the MAX3 test to ex, asy method is used to approximate the p-value ##of MAX3
## "MAX3 test"
##$method
## "asy"
## the method chosen is asy
##$statistics
## 5.117125
## The statistics of the test is 5.117125
##$Pvalue
## 8.562347e-07
## The associated p-value of the test using asy method is 8.562347e-07
##$conclusion
## "alternative hypothesis: association exists under significant level 0.05"
## Under the significant level 0.05, the alternative hypothesis that the 

##association exists is accepted according to the given p-values.

}
\seealso{\code{\link{CATT}}
         \code{\link{caco}}}

