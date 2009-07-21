\name{ABT}
\Rdversion{1.1}
\alias{ABT}
\title{
The allelic based test
}
\description{
This function conducts the allelic based test to a 2 by 3 case-control contingency table and reports the test statistics, associated p-value and conclusion for the hypothesis test respectively.
}
\usage{
ABT(data)
}

\arguments{
  \item{data}{

data is a 2 by 3 case-control contingency table. The first and second rows represent the case group and control group respectively. The first, second and third columns represent the genotypes of a susceptibility diallelic marker containing 0, 1 and 2 risk allele respectively. Thus, the numbers in this table represent the genotype counts belonging to the corresponding genotypes and case-control status. 
}
}
\value{
\item{statistics}{

the statistics of the allelic based test}
\item{Pvalue}{

the associated p-value of the allelic based test}
\item{conclusion}{

the conclusion of the hypothesis test}
}
\references{
Sasieni, PD(1997). From genotypes to genes: doubling the sample size. Biometrics 53, 1253-1261.

Zang Y, Fung WK and Zheng G(2009). Simple Algorithms to calculate the asymptotic null distribution for MAX3 and genetic model selection in case-control genetic association studies in R. Journal of Statistical software (revised).  
}
\examples{
library(Rassoc)
data(caco)
ex=matrix(caco[1,],nrow=2,byrow=TRUE)
## ex is an example of a 2 by 3 case-control contingency table. 
ABT(ex)
## Conduct the allelic based test to dataset ex. 
##"the allelic based test"
##$statistics
##-5.490279
## The statistics of the test is -5.490279.
##$Pvalue
##4.012992e-08
## The associated p-value of the test is 4.012992e-08.
##$conclusion
##"alternative hypothesis: association exists under significant level 0.05"
## Under the significant level 0.05, the alternative hypothesis that the


##association exists is accepted according to the given p-values.
}
\seealso{\code{\link{caco}}}

