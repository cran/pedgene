## Purpose: Test script for pedgene package
## Testing a simple 10-variant gene on both the X chrom and autosome
## Authors: Dan Schaid and Jason Sinnwell
## Created: 19-AUG-2013
## Updated: 11/7/2013

## alternatively:
require(kinship2)
require(pedgene)

data(example.ped)
data(example.geno)
data(example.map)


if(0) {
## delete above line to look at pedigrees
## quick look at the 3 pedigrees
  pedall <- with(example.ped, pedigree(famid=ped, id=person, dadid=father,
                momid=mother, sex=sex, affected=ifelse(is.na(trait), 0, trait)))
  
  plot(pedall[1])
  plot(pedall[2])
  plot(pedall[3])
}


## simple tests of two genes (10 variants each)
## the genes are same variants, just on chroms 1 and X
pg.m2 <- pedgene(example.ped, example.geno, example.map, male.dose=2)

pg.m1 <- pedgene(example.ped, example.geno, example.map, male.dose=1)


## saved objects are of class pedgene, with items call (function call)
## and pgdf, a data.frame with a row for each gene
class(pg.m2)
names(pg.m2)


print(pg.m2, digits=4)


print(pg.m1, digits=4)

summary(pg.m2)

## Testing first gene with dose=2-dose
geno.recode <- cbind(example.geno[,1:2], 2-example.geno[,grep("AA", names(example.geno))])
pg.recode <- pedgene(example.ped, geno.recode, male.dose=2)

## note when map not given, assumes all 1 gene, and assigns "unknown" gene/chrom
pg.recode
