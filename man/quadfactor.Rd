\name{quadfactor}
\alias{quadfactor}
\title{
  Internal utility to calculate a constant quadratic factor for
  gene-level statistics variances, over all pedigrees
}
\description{
  Internal utility to calculate a constant quadratic factor for
  gene-level statistics
  variances over all pedigrees, for either autosomes or X chromosome
}
\usage{
quadfactor(kinmat, chrom, resid, sex, male.dose)
}
\arguments{
  \item{kinmat}{A kinship matrix for one or more pedigrees}
  \item{chrom}{character string for chromosome number, if "X", the
    method accounts for sex code}
  \item{resid}{the residual based on the trait minus the group mean or
    adjusted trait value for each subject}
  \item{sex}{See pedgene}
  \item{male.dose}{See pedgene}
}
\value{
  Constant quadratic factor for gene-level statistics variances, for
  autosomes and X chromosome 
}
\author{
Daniel J. Schaid, Mayo Clinic (schaid@mayo.edu).
}
\seealso{
pedgene, kinship
}

