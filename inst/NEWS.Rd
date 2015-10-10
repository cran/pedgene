\name{NEWS}
\title{NEWS file for the pedgene package}
\section{Release Notes for version 2.9}{
  \itemize{
    \item Use trait.adjusted as mean of trait for those with genotype data
  }
}
\section{Release Notes for version 2.6}{
  \itemize{
    \item Make sure character ids stay as character ids
    \item Bug fix for no relateds
    \item For X chrom, add check for monomorphic to include within males
    and females only
  }
}
\section{Release Notes for version 1.4}{
  \itemize{
    \item Add argument acc.davies for pedgene. Some kernel tests gave NA
    p-values with default daves accuracy of 1e-6. Performs well with
    default of 1e-5, but give option to user in case NAs occur.
    \item return n-variants for the gene tests, now a column in the returned
    object from pedgene
    \item Allow the method to complete if any gene has no markers even
    after all are removed for having zero variance (same genotype)
    \item When a gene has only one marker, the kernel test reduces to
    burden test. Return same test stat and p-value for both. Thanks to
    Seth Ament.
    \item Fix bug where failed to subset genotype matrix when either NA
    in trait or all NA genotypes for a given person, after matching geno
    and ped subject by ped-person ids. Thanks again to Seth Ament. 
    \item Change maintainer to Jason, rather than Dan.
    \item Add n.noninform to pedgene object, the number of
    non-informative markers that were removed
    
  }
}
\section{Release Notes for version 1.0}{
  \itemize{
    \item Require ped and person ID columns in genotype data. To match to
    more-complete subjects in pedigree data by ped and person ID.
    \item Added argument checkpeds, allows user to turn off the step
    that performs pedigree checks.
    \item Kinship function (from kinship2) is now called without making
    a pedigree object first. This step allows disconnected pedigrees, and
    identifiers can be either numberic or character.
    \item Perform Davies' p-value for Kernel statistic.
    \item Added package vignette, a simple manual
  }
}    
\section{Release Notes for version 0.5}{
  \itemize{
    \item wrapper called pedgene to call pedgene.stats on one gene at a
    time while re-using kinship matrix for all pedigrees    
    \item require ped and geno data.frames to have same number of rows
    \item add test script: test.pedgene.R
    \item if map given (with chrom and gene columns), require map to have
    as many rows as geno has columns
    \item define S3 pedgene object, with simple print and plot methods
  }
}
\section{Release Notes for version 0.1}{
\itemize{
  \item Original release, mostly code from Dan's scripts
  \item pedgene.init pre-calculates some stats, passed to pedgene.stats
  for one single gene
}
}
