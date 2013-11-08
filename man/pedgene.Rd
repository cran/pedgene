\name{pedgene}
\alias{pedgene}
\alias{pedgene.stats}
\alias{print.pedgene}
\alias{summary.pedgene}
\title{
 Compute Kernel and Burden Statistics for Pedigree Data (possibly
 with unrelated subjects) 
}
\description{
  Compute linear kernel and burden statistics for gene-level analysis of
  data that includes pedigree-related subjects, and possibly unrelated
  subjects.
}
\usage{
pedgene(ped, geno, map=NULL, male.dose=2, weights=NULL, checkpeds=TRUE)
}
\arguments{
  \item{ped}{A data.frame with variables that define the pedigree structure
    (typical format used by LINKAGE and PLINK), trait (phenotype), and
    optionally a covariate-adjusted trait (for covariate-adjusted gene
    level statistics). The variables in the data.frame must be named as
    follows:
    \itemize{
      \item{ped: pedigree ID, character or numeric allowed}
      \item{person: person ID, a unique ID within each pedigree, numeric
      or character allowed}
      \item{father: father ID, 0 if no father}
      \item{mother: mother ID, 0 if no mother}
      \item{sex: coded as 1 for male, 2 for female}
      \item{trait: phenotype, either case-control status coded as 1
	for affected and 0 for unaffected, or a continuous
	value. Subjects with missing (NA) will be removed from the analysis}
      \item{trait.adjusted: an optional variable for covariate-adjusted
        trait. If trait.adjusted is present in the data.frame, then
	gene-level tests are adjusted for covariates using
        residuals = (trait - trait.adjusted). Otherwise, gene-level tests
        are not adjusted for covariates, in which case residuals = trait - mean(trait)}
    }
  }
  \item{geno}{
    Data.frame or matrix with genotypes for subjects (rows) at each variant position
    (columns). The first two columns are required to be named \code{ped} and
    \code{person}, which are used to match subjects to their data in the
    \code{ped} data.frame. The genotypes are coded as 0, 1, 2 for autosomal
    markers (typically a count of the number of the less-frequent
    allele). For X-chromosome markers, females are coded 0, 1, 2, and
    males coded 0, 1. Missing genotypes (NA) are allowed.
  }
  \item{map}{
    Optional data.frame with columns "chrom" and  "gene", one row per
    variant column in geno. The gene name can be any identifier for the
    gene. The chromosome can be either numeric or character, where the
    calculations will differ between autosomes vs X chromosome (allow
    "X"/"x"/23, converted to "X" in results)
  }
  \item{male.dose}{When analyzing the X-chromosome, male.dose defines
    how male genotypes should be analyzed. male.dose can be between 0
    and 2, but is typically either 1 or 2. Ozbek and Clayton show that
    male.dose = 2 is powerful in the presence of X-chromosome dosage
    compensation in females.
  }
  \item{weights}{optional user-specified weights, a vector of weights
    for each variant column of geno (M). If none given, Madsen-Browning
    weights are applied, where these weights are 1/sqrt(maf*(1-maf))}
  \item{checkpeds}{logical, if FALSE, the method will skip the pedigree
    checking step, which can be intensive for large studies}
}
\details{
  The pedgene function is a wrapper function to call pedgene.stats on
  one gene at a time. The pedgene.stats function calculates gene-level
  tests for associations with a trait among subjects, accounting for relationships
  among subjects based on known pedigree relationships. This is achieved
  by the kinship function in the kinship2 package. The kernel
  association statistic uses a weighted linear kernel, with default
  weights based on those suggested by Madsen
  and Browning (weight = 1 / sqrt(maf*(1-maf)), where maf is the minor
  allele frequency for a variant). The burden statistic is based on a
  weighted sum of variants, also using the Madsen-Browning weights.
  Variant positions that have zero variance are removed from the
  analysis because they do not contribute information.
}
\value{
  An object of the pedgene S3 class, with the following elements:
  \item{call: }{function call}
  \item{pgdf: }{data.frame with gene name, chromosome, kernel and burden
    test statistics and p-values. Kernel p-values are based on Davies'
    method (Davies 1980, reference below), and burden p-values are
    based on the normal distribution}
}
\references{
Schaid DJ, McDonnell SK., Sinnwell JP, Thibodeau SN. (2013)
Multiple Genetic Variant Association Testing by Collapsing and Kernel
Methods With Pedigree or Population Structured Data. Genetic Epidemiology, 37(5):409-18.

Davies R.B. (1980) Algorithm AS 155: The Distribution of a Linear Combination
of chi-2 Random Variables, Journal of the Royal Statistical
Society. Series C (Applied Statistics), 29(3):323-33
}
\author{
Daniel J. Schaid, Jason P. Sinnwell, Mayo Clinic (schaid@mayo.edu).
}
\seealso{
pedigreeChecks, example.ped
}
\examples{
# example data with the same 10 variants for an autosome and X chromosome
# pedigree data on 39 subjects including 3 families and unrelateds
data(example.ped)
data(example.geno)
data(example.map)

# gene tests (chroms 1 and X) with male.dose=2
pg.m2 <- pedgene(example.ped, example.geno, example.map, male.dose=2)
# same genes, with male.dose=1
pg.m1 <- pedgene(example.ped, example.geno, example.map, male.dose=1)

## print and summary methods
print(pg.m2, digits=3)
summary(pg.m1, digits=3)
}
\keyword{kinship}

