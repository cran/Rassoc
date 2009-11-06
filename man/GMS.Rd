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
\item{statistic}{

the statistic of the GMS test}
\item{p.value}{

the associated p-value of the GMS test}
}
\references{
Joo, J, Kwak, M and Zheng, G (2009). Improving power for testing genetic association in case-control studies by reducing alternative space. Biometrics (in press).

Zang Y, Fung WK and Zheng G(2009). Simple Algorithms to calculate the asymptotic null distribution for MAX3 and genetic model selection in case-control genetic association studies in R. Journal of Statistical software (revised).

Zheng, G and Ng HKT (2008). Genetic model selection in two-phase analysis for
case-control association studies. Biostatistics 9, 391-399.
}
\examples{
library(Rassoc)
ca=c(139,249,112) 
co=c(136,244,120) 
a=rbind(ca,co) 

GMS(a,"boot",100000) 
## The GMS test using the boot method
## data:  a 
## statistic = 0.4894, p-value = 0.6658

GMS(a,"bvn",100000)
## The GMS test using the bvn method
## data:  a 
## statistic = 0.4894, p-value = 0.6598

GMS(a,"asy",100000)
## The GMS test using the asy method
## data:  a 
## statistic = 0.4894, p-value = 0.6621
}
\seealso{\code{\link{CATT}}}

