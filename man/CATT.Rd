\name{CATT}
\Rdversion{1.1}
\alias{CATT}
\title{
The Cochran-Armitage trend test
}
\description{
This function conducts the Cochran-Armitage trend test to a 2 by 3 case-control contingency table and reports the test statistics, associated p-value and conclusion for the hypothesis test respectively.
}
\usage{
CATT(data, x)
}
\arguments{
  \item{data}{
data is a 2 by 3 case-control contingency table. The first and second rows represent the case group and control group respectively. The first, second and third columns represent the genotypes of a susceptibility diallelic marker containing 0, 1 and 2 risk allele respectively. Thus, the numbers in this table represent the genotype counts belonging to the corresponding genotypes and case-control status. 
}
  \item{x}{
x is the score of the Cochran-Armitage trend test. It can be any real number between 0 and 1. Specifically, x=0, 0.5 and 1 are optimal for recessive, additive/multiplicative and dominant genetic models respectively. 
}
}
\value{
\item{statistics}{

the statistics of the Cochran-Armitage trend test}
\item{Pvalue}{

the associated p-value of the Cochran-Armitage trend test}
\item{conclusion}{

the conclusion of the hypothesis test}
}
\references{
Sasieni, PD(1997). From genotypes to genes: doubling the sample size. Biometrics 53, 1253-1261.

Zang Y, Fung WK and Zheng G(2009). Simple Algorithms to calculate the asymptotic null distribution for MAX3 and genetic model selection in case-control genetic association studies in R. Journal of Statistical software (revised).

Zheng, G, Freidlin, B, Li, Z and Gastwirth, JL (2003). Choice of scores in trend tests for case-control studies of candidate-gene associations. Biometrical Journal 45, 335-348. 
}
\examples{
library(Rassoc)
data(caco)
ex=matrix(caco[1,],nrow=2,byrow=TRUE)
## ex is an example of a 2 by 3 case-control contingency table. 
CATT(ex,0.5)
## Conduct the Cochran-Armitage trend test to dataset ex. x is chosen as 0.5, ##which is optimal under the additive/multiplicative genetic model.
## "the Cochran-Armitage trend test"
##$statistics
## -5.117125
## The statistics of the test is -5.117125.
##$Pvalue
## 3.102276e-07
## The associated p-value of the test is 3.102276e-07.
##$conclusion
## "alternative hypothesis: association exists under significant level 0.05"
## Under the significant level 0.05, the alternative hypothesis that the 

##association exists is accepted according to the given p-values.

}
\seealso{\code{\link{caco}}}

