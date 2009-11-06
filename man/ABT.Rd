\name{ABT}
\Rdversion{1.1}
\alias{ABT}
\title{
The allelic based test
}
\description{
This function conducts the allelic based test to a 2 by 3 case-control contingency table and reports the test statistics and associated p-value.}
\usage{
ABT(data)
}

\arguments{
  \item{data}{

data is a 2 by 3 case-control contingency table. The first and second rows represent the case group and control group respectively. The first, second and third columns represent the genotypes of a susceptibility diallelic marker containing 0, 1 and 2 risk allele respectively. Thus, the numbers in this table represent the genotype counts belonging to the corresponding genotypes and case-control status. 
}
}
\value{
\item{statistic}{

the statistic of the allelic based test}
\item{p.value}{

the associated p-value of the allelic based test}
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
## The allelic based test
## data:  ex 
## statistic = -5.4903, p-value = 4.013e-08
## The statistic of the test is -5.4903.
## The associated p-value of the test is 4.013e-08.
}
\seealso{\code{\link{caco}}}

