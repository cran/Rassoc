\name{GMS}
\Rdversion{1.1}
\alias{GMS}
\title{
The genetic model selection test
}
\description{
This function conducts the genetic model selection (GMS) test to a 2 by 3 case-control contingency table using either empirical or asymptotic methods. GMS is a two-phase test. In phase 1, the Hardy-weinberg disequilibrium trend test (HWDTT) is used to detect the underlying genetic model. In phase 2, an optimal Cochran-Armitage trend test corresponding to the selected genetic model is used for testing association.  
}
\usage{
GMS(data, method, m)
}
\arguments{
  \item{data}{
data is a 2 by 3 case-control contingency table. The first and second rows represent the case group and control group respectively. The first, second and third columns represent the genotypes of a susceptibility diallelic marker containing 0, 1 and 2 risk allele respectively. Thus, the numbers in this table represent the genotype counts belonging to the corresponding genotypes and case-control status. 
}
  \item{method}{
method is the approach used for calculating the p-value of GMS statistics. boot and bvn are empirical methods and asy is asymptotic method.  
}
  \item{m}{
m is the replication times used for calculating the empirical p-values (boot and bvn). It can be any positive integer (1, for example) while choosing asymptotic method (asy). 
}
}
\details{
boot method re-samples case-control data under null hypothesis and applies GMS to each replicate. bvn directly generates statistics of GMS rather than case-control data under null hypothesis. asy calculates the p-value from the asymptotic distribution of GMS under null hypothesis. 
}
\value{
\item{method}{

the method chosen to calculate the p-value of GMS}
\item{statistics}{

the statistics of the GMS test}
\item{Pvalue}{

the associated p-value of the GMS test}
\item{conclusion}{

the conclusion of the hypothesis test}

}
\references{
Joo, J, Kwak, M and Zheng, G (2009). Improving power for testing genetic association in case-control studies by reducing alternative space. Biometrics (in press).

Zang Y, Fung WK and Zheng G(2009). Simple Algorithms to calculate the asymptotic null distribution for MAX3 and genetic model selection in case-control genetic association studies in R. Journal of Statistical software (revised).

Zheng, G and Ng HKT (2008). Genetic model selection in two-phase analysis for
case-control association studies. Biostatistics 9, 391-399.
}
\examples{
library(Rassoc)
data(caco)
ex=matrix(caco[1,],nrow=2,byrow=TRUE)
## ex is an example of a 2 by 3 case-control contingency table. 

## Using boot method 
GMS(ex,"boot",10000)
## Conduct the GMS test to ex, boot method is used based on 1000000   
##replicates to approximate the p-value of GMS
## "genetic model selection test"
##$method
## "boot"
## the method chosen is boot
##$statistics
## 5.117125
## The statistics of the test is 5.117125
##$Pvalue
## 2e-06
## The associated p-value of the test using boot method is 2e-06
##$conclusion
## "alternative hypothesis: association exists under significant level 0.05"
## Under the significant level 0.05, the alternative hypothesis that the 

##association exists is accepted according to the given p-values.

## Using bvn method 
GMS(ex,"bvn",10000)
## Conduct the GMS test to ex, bvn method is used based on 1000000   
## replicates to approximate the p-value of GMS
## "genetic model selection test"
##$method
## "bvn"
## the method chosen is bvn
##$statistics
## 5.117125
## The statistics of the test is 5.117125
##$Pvalue
## 1e-06
## The associated p-value of the test using bvn method is 1e-06
##$conclusion
## "alternative hypothesis: association exists under significant level 0.05"
## Under the significant level 0.05, the alternative hypothesis that the 

##association exists is accepted according to the given p-values.

## Using asy method 
GMS(ex,"asy",1)
## Conduct the GMS test to ex, asy method is used to approximate the p-value ##of GMS
## "genetic model selection test"
##$method
## "asy"
## the method chosen is asy
##$statistics
## 5.117125
## The statistics of the test is 5.117125
##$Pvalue
## 8.61515e-07
## The associated p-value of the test using asy method is 8.61515e-07
##$conclusion
## "alternative hypothesis: association exists under significant level 0.05"
## Under the significant level 0.05, the alternative hypothesis that the 

##association exists is accepted according to the given p-values. 
}
\seealso{\code{\link{CATT}}
         \code{\link{caco}}}

