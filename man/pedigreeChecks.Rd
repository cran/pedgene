\name{pedigreeChecks}
\alias{pedigreeChecks}
\title{
Validity checks on pedigree data
}
\description{
Checks for valid IDs, sex codes for data in a single pedigree
}
\usage{
pedigreeChecks(pedigree, male.code = 1, female.code = 2)
}
\arguments{
  \item{pedigree}{data frame with variables named person, father, mother, sex.}
  \item{male.code}{sex code for males}
  \item{female.code}{sex code for females}
}
\details{
 A series of basic pedigree checks
}
\value{
valid = TRUE or FALSE for validity of pedigree data
}
\author{
Daniel J. Schaid (schaid@mayo.edu).
}

