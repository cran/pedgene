
## Functions: pedgene 
## wrapper for computing retrospective likelihood stat on pedigrees for rare
## variants over multiple genes
## Authors: Jason Sinnwell and Dan Schaid

pedgene <- function(ped, geno, map=NULL, male.dose=2, weights=NULL, checkpeds=TRUE,
                    acc.davies=1e-5) {

##Arguments:
##  
##  Required:
##  ped: data.frame with columns needed to create the pedigree
##  geno: data.frame with ped, person ids in the first two columns, and numeric columns
##       with minor allele count (0/1/2) for markers (columns) and subjects (rows) in
##       the ped object
##  
##  Optional:
##  map: data.frame with columns chrom, position, and gene.
##        gene can be gene symbol or geneid (entrez or ensemble);
##        it is not matched with other annotation, just used for marker groups
##        If not passed, assume all variants from the same gene 
##  male.dose: When doing X-chrom, define how male genotypes should be
##        analyzed. male.dose can be between 0 and 2, but usually 1 or 2
##  weights: allow user-specified weights. If Null, use Madsen-Browning weights

   
## Steps
## 1) verify ped columns, map and geno match dimensions
## 2) Create kinship matrices for autosomes and X for all subjects
## 3) run pedgene.stats on each gene given in map

  verbose=FALSE
  
  ## save the call
  call <- match.call() 
  ## save options before setting stringsAsFactors for just this function
  saveOpt <- options()
  options(stringsAsFactors=FALSE)
  
  ## require kinship function to be recent
  kin2v <- sessionInfo()$otherPkgs$kinship2$Version
  if(is.null(kin2v)) {
    kin2v <- sessionInfo()$loadedOnly$kinship2$Version
  }
  if(as.numeric(substring(kin2v, 1, nchar(kin2v)-2)) < 1.5) {
    stop("kinship2 needs to be version 1.5.3 or later\n")
  }
  
  ## if no map, create one, assuming all one gene
  if(is.null(map)) {
    map <- data.frame(chrom=rep("unknown", ncol(geno)-2),
                      gene=rep("unknown", ncol(geno)-2))
  }

  ## unify names of ped and map to lowercase
  names(map) <- casefold(names(map))
  names(ped) <- casefold(names(ped))

  ## verify map data.frame ###
  ## i. check column names
  if(any(!(c("chrom", "gene") %in% names(map)))) {
    stop("map requires columns for chrom and gene")
  }
  ## ii. recode chrom 23 to X
  map$chrom[map$chrom==23] <- "X"

  ## verify geno matrix, and that it matches map
  if(any(!(c("ped", "person") %in% names(geno)))) {
    stop("geno requires columns 'ped' and 'person' ids")
  }
  ## get indices of ped/person of geno to match ped, then strip off those columns
  keepped <- match(paste(geno$ped, geno$person, sep="-"),
                   paste(ped$ped, ped$person, sep="-"))
  geno <- geno[,!(names(geno) %in% c("ped", "person"))]
  
  if(nrow(map) != (ncol(geno))) {
    stop(paste("map rows (", nrow(map), ") and geno columns (", ncol(geno),
               ") do not match \n",sep=""))
  }
  

  ## Check that geno for males on X should only have 0 and 1 dosages
  xidx <- which(map$chrom=="X" | map$chrom=="x")
  if(length(xidx)) {
    xdosemale <- geno[ped$sex[keepped]==1,xidx, drop=TRUE]
    if(sum(xdosemale>1, na.rm=TRUE)) {
      stop("All male dose on X chromosome should be <= 1")
    }
  }
  
  ## verify ped data.frame has expected column names
  if(any(!(c("ped", "person", "father", "mother", "sex", "trait")
           %in% names(ped)))) {
    stop("Error: ped requires columns: ped, person, father, mother, sex, trait")
  }
   
  ## add trait.adjusted if not already there
  if(!("trait.adjusted" %in% names(ped))) {
    ped$trait.adjusted <- mean(ped$trait, na.rm=TRUE)      
  }
  
  ## verify weights, match ncol(geno)
  if(!is.null(weights)) {
    ## by default, do Madsen-Browning weights, currently within ped.gene.stats
    ## otherwise, these are user-specified, check length
    if(length(weights) != ncol(geno)) {
       stop(paste("Error: should have weights(", length(weights),
                  ") for every variant position(", ncol(geno), ")", sep=""))
    }
  }

  ## perform simple pedigree checks
  if(checkpeds) {
    uped <- unique(ped$ped)
    nped <- length(uped)
    
    for(i in 1:nped) {
      
      iped <- uped[i]
      
      temp.ped <- ped[ped$ped == iped,, drop=FALSE]
      
      if(nrow(temp.ped) > 1) {      
        ## simple checks on pedigree
        pedigreeChecks(temp.ped, male.code=1, female.code=2)
      }
    }
  }
   ## additional checks could be done on peds when creating pedlist object,
  ## which could be used to create kinmat. We rather created it directly from ped
  # pedall <- with(ped, kinship2::pedigree(id=person, dadid=father, momid=mother,
  #                          sex=sex, famid=ped, missid=missid))
  #  kinmat <- kinship2::kinship(pedall, chrtype="auto")
    

  ## create kinship matrix, also for X if any genes on X chrom
  ## subset to only those with geno rows
  kinmat <- Matrix(with(ped, kinship(ped, id=paste(ped,person,sep="-"),
                      dadid=ifelse(father>0,paste(ped,father,sep="-") , as.character(father)),
                      momid=ifelse(mother>0, paste(ped,mother,sep="-"), as.character(mother)),
                      sex=sex, chrtype="autosome")))
  kinmat <- kinmat[keepped, keepped]
 
  if(any(map$chrom=="X")) {
    kinmatX <- Matrix(with(ped, kinship(ped, id=paste(ped,person,sep="-"),
                      dadid=ifelse(father>0,paste(ped,father,sep="-") , as.character(father)),
                      momid=ifelse(mother>0, paste(ped,mother,sep="-"), as.character(mother)),
                      sex=sex, chrtype="X")))
    kinmatX <- kinmatX[keepped, keepped]
  } else {
    kinmatX <- NULL
  }
  ped <- ped[keepped,]
  
  ## subset pedgeno kinmat, kinmatX to only subject who have genotype data
  missidx <- is.na(ped$trait) | apply(is.na(geno), 1, all) 
  if("trait.adjusted" %in% names(ped)) missidx <- missidx | is.na(ped$trait.adjusted)
  if(sum(missidx)>0) {
    ped <- ped[!missidx,]
    kinmat <- kinmat[!missidx, !missidx]
    kinmatX <- kinmatX[!missidx, !missidx]
    geno <- geno[!missidx,]
  }
  
  gvec <- chromvec <- nvariant <- kstat <- kpvaldav <- bstat <- bpval <- NULL
  
  for(g in unique(map$gene)) {
    if(verbose) {
      cat("test on gene ", g, "\n")
    }
    gidx <- which(map$gene==g)
    ## JPS add drop=FALSE for 1-marker genes 11/13/13
    genosub <- geno[,gidx,drop=FALSE]

    resid <- ped$trait - ped$trait.adjusted
    sex <- ped$sex
    chrom <- map$chrom[gidx[1]]
    
    c.factor <- quadfactor(
             if(chrom=="X") kinmatX else kinmat,
             chrom, resid, sex, male.dose)

    pgstat <- pedgene.stats(genosub, as.vector(c.factor), map$chrom[gidx[1]],
                    male.dose, sex, resid, weights=weights[gidx], acc.davies=acc.davies)
    if(pgstat$nvariant==0) {
      cat("gene: '", g, "' has no markers after removing markers with all same genotype\n")
    }
    gvec <- c(gvec, g)
    chromvec <- c(chromvec, chrom)
    nvariant <- c(nvariant,pgstat$nvariant)
    kstat <- c(kstat, pgstat$stat.kernel)
    kpvaldav <- c(kpvaldav, pgstat$pval.kernel.davies)
    bstat <- c(bstat, pgstat$stat.burden)
    bpval <- c(bpval, pgstat$pval.burden)   
  }
  
  pgdf <- data.frame(gene=gvec, chrom=chromvec, nvariant=nvariant, stat.kernel=kstat,
                     pval.kernel.davies=kpvaldav,
                     stat.burden=bstat, pval.burden=bpval)
  
  # re-set options
  options(saveOpt)
  
  pglist <- list(pgdf=pgdf, call=call)
  class(pglist) <- "pedgene"
  return(pglist)

}

## print and summary methods for pedgene S3 class

print.pedgene <- function(x, ...) {
 
  print.data.frame(x$pgdf, ...)

  invisible()
}
summary.pedgene <- function(object, ...) {
  
  cat("\nSummary for pedgene object: \n\n")
  cat("Call:\n")
  print(object$call)
  cat("\n\n")

  ## invoke print method
  print(object,  ...)
  
  invisible()
}
