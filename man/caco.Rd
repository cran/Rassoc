\name{caco}
\Rdversion{1.1}
\alias{caco}
\docType{data}
\title{
The case-control dataset of 17 SNPs
}
\description{
caco is a matrix containing the case-control samples of 17 SNPs reported from four genome-wide association studies for age-related macular degeneration (AMD) (Klein et al., 2005), cancer studies (Hunter et al., 2007 and Yeager et al., 2007) and a hypertension study (WTCCC, 2007).
}
\usage{data(caco)}
\details{
caco is a 17 by 6 matrix. Each row represents a different SNP. The 1st, 2nd and 3rd columns represent the genotypes of case group containing 2, 1 and 0 risk alleles; the 4th, 5th and 6th columns represent the genotypes of control group containing 2, 1 and 0 risk alleles. Thus, each number in the caco matrix represents the genotype counts belonging to the corresponding SNP, genotype and case-control status. 
}
\source{
Hunter DJ, Kraft P, Jacobs KB, Cox DG, Yeager N, Hankinson SE, Wacholder S, Wang Z, Welch R, Hutchinson A, et al (2007). A genome-wide association study identifies alleles in FGFR2 associated with risk of sporadic postmenopausal breast cancer. Nature Genetics 39, 870-874.

Klein RJ, Zeiss C, Chew EY, Tsai J-Y, Sackler RS, Haynes C, Henning AK, SanGiovanni JP, Mane SM, Mayne ST, Bracken MB, Ferris FL et al (2005). Complement factor H polymorphism in aged-related macular degeneration. Science 308, 385-389.

The Wellcome Trust Case Control Consortium (WTCCC) (2007). Genome-wide association study of 14,000 cases of seven common diseases and 3,000 shared controls. Nature 447, 661-683.

Yeager M, Orr N, Hayes RB, Jacobs KB, Kraft P, Wacholder S, Minichiello MJ, Fearnhead P, Yu K, Chatterjee N et al (2007). Genome-wide association study of prostate cancer identifies a second risk locus at 8q24. Nature Genetics 39, 645-649.
}
\examples{
data(caco)
ex=matrix(caco[1,],nrow=2,byrow=TRUE)
## ex is a 2 by 3 case-control contingency table of SNP rs380390.
ABT(ex)
## Conduct the allelic based test to dataset ex. 
## "the allelic based test"
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

