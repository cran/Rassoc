\name{MERT}
\Rdversion{1.1}
\alias{MERT}
\title{
The maximin efficiency robust test
}
\description{
This function conducts the maximin efficiency robust test to a 2 by 3 case-control contingency table and reports the test statistics, associated p-value and conclusion for the hypothesis test respectively.
}
\usage{
MERT(data)
}
\arguments{
  \item{data}{
data is a 2 by 3 case-control contingency table. The first and second rows represent the case group and control group respectively. The first, second and third columns represent the genotypes of a susceptibility diallelic marker containing 0, 1 and 2 risk allele respectively. Thus, the numbers in this table represent the genotype counts belonging to the corresponding genotypes and case-control status. 
}
}
\value{
\item{statistics}{

the statistics of the maximin efficiency robust test}
\item{Pvalue}{

the associated p-value of the maximin efficiency robust test}
\item{conclusion}{

the conclusion of the hypothesis test}
}
\references{
Freidlin, B, Zheng, G, Li, Z and Gastwirth, JL (2002). Trend tests for case-control studies of genetic markers: power, sample size and robustness. Human Heredity 53, 146-152.

Gastwirth JL (1966). On robust procedures. Journal of American Statistical Association 61, 929-948. 

Gastwirth JL (1985). The use of maximin efficiency robust tests in combining contingency tables and survival analysis. Journal of American Statistical Association 80, 380-384. 

Zang Y, Fung WK and Zheng G(2009). Simple Algorithms to calculate the asymptotic null distribution for MAX3 and genetic model selection in case-control genetic association studies in R. Journal of Statistical software (revised).  

}
\examples{
library(Rassoc)
data(caco)
ex=matrix(caco[1,],nrow=2,byrow=TRUE)
## ex is an example of a 2 by 3 case-control contingency table. 
MERT(ex)
## Conduct the maximin efficiency robust test to dataset ex. 
## "the maximin efficiency robust test"
##$statistics
## -5.073485
## The statistics of the test is -5.073485.
##$Pvalue
## 3.905944e-07
## The associated p-value of the test is 3.905944e-07.
##$conclusion
## "alternative hypothesis: association exists under significant level 0.05"
## Under the significant level 0.05, the alternative hypothesis that the 

##association exists is accepted according to the given p-values.

}
\seealso{\code{\link{caco}}}

