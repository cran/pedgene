\name{example.geno}
\alias{example.geno}
\alias{example.ped}
\alias{example.map}
\alias{example.relation}
\docType{data}
\title{Example datasets for pedgene }
\description{
  example.geno: a data frame with minor allele count for subjects (rows)
  at variant positions (columns); example.ped: pedigree and trait data for subjects in
  example.geno; example.map: gene and chromosome for variant positions in example.geno;  
  example.relation: special (twin) relationships for individuals in example.ped
}
\usage{
data(example.geno)
data(example.ped)
data(example.map)
data(example.relation)
}
\format{
  \code{example.geno} -data frame with minor allele count for 20 variant
  positions:
  \describe{
    \item{\code{ped }}{pedigree ID, character or numeric}
    \item{\code{person }}{person ID, used with ped to match subjects to
	their row in example.ped}
    \item{\code{AA.1-AA.10,AX.1-AX.10 }}{genotype columns at 10 positions for
      each of 2 simulated genes}
  }
  \code{example.ped} -data frame with pedigree structure and trait
  values in the following columns:
  \describe{
    \item{\code{ped }}{ pedigree ID, character or numeric}
    \item{\code{person }}{ person ID, a unique ID within each pedigree}
    \item{\code{father }}{father ID, 0 if no father}
    \item{\code{mother }}{mother ID, 0 if no mother}
    \item{\code{sex }}{coded as 1 for male, 2 for female}
    \item{\code{trait }}{phenotype, either case-control status coded as 1
      for affected and 0 for unaffected, or a continuous
      value. Subjects with missing (NA) will be removed from the analysis}
    \item{\code{trait.adjusted }}{an optional variable for covariate-adjusted
        trait. If trait.adjusted is present in the data.frame, then
	gene-level tests are adjusted for covariates using
        residuals = (trait - trait.adjusted). Otherwise, gene-level tests
        are not adjusted for covariates, in which case residuals = trait - mean(trait)}
    }
 \code{example.map } - data frame with columns for gene name and chromosome:
 \describe{
   \item{\code{chrom }}{chromosome code (1-23,X allowed) where the gene is located}
   \item{\code{gene }}{gene identifier}
 }
 \code{example.relation } - matrix specifying special relationships between pairs 
 	of individuals in the following columns:
    \describe{
      \item{\code{ped }}{pedigree ID, character or numeric}
      \item{\code{id1 }}{person1 ID, used with ped to match subjects to
	their row in example.ped}
      \item{\code{id2 }}{person2 ID, used with ped to match subjects to
	their row in example.ped}
      \item{\code{code }}{any of the following: 1=Monozygotic twin, 2=Dizygotic twin, 
	  		3=Twin of unknown zygosity, and 4=Spouse}
    }
}
\source{
 Simulated data for testing purposes
}
\examples{
data(example.geno)
data(example.ped)
data(example.map)
data(example.relation)
}
\keyword{datasets}
