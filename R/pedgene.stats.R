
## Title: pedgene.stats.R
## Purpose: worker function to calculate ped-based kernel and burden tests for rare variants
##          on one gene at a time. Wrapper function over multiple genes is pedgene()
## Author: Jason Sinnwell
## Created: 8/2013
## Updated: 10/30/2013


pedgene.stats <- function(geno, c.factor, chrom, male.dose, sex,
                          resid, weights=NULL, accuracy.davies=1e-6) {

    x.chrom <- as.character(chrom)=="X"
    
    # remove monomorphic markers
    v <- apply(geno, 2, var, na.rm=TRUE)
    geno <- geno[, v > 0]

    geno <- as.matrix(geno)
    if(ncol(geno) == 0) {
        stop("No data in geno after remove monomorphic markers")
    }
    
    # Account for missing genotypes. Could exclude subjects with
    # any missing, but  here just replace missing with mean values.
    # Would be best for user to decide how to handle missing prior
    # to calling this function.

    col.miss <- apply(is.na(geno), 2, any)
   
    if(any(col.miss == TRUE))
      {
        for(j in 1:ncol(geno))
          {
            if(col.miss[j])
              {
                if(x.chrom)
                  {
                    mn.male <- mean(geno[sex==1,j], na.rm=TRUE)
                    geno[is.na(geno[,j])&sex==1, j] <- mn.male
                    mn.female <- mean(geno[sex==2,j], na.rm=TRUE)
                    geno[is.na(geno[,j])&sex==2, j] <- mn.female          
                  } else
                {
                  mn <- mean(geno[,j], na.rm=TRUE)
                  geno[is.na(geno[,j]), j] <- mn
                }
              }
          }
      }
   
    # Madsen-Browning weights
    if(x.chrom)
      {
        count.male   <- apply(geno[sex==1,], 2, sum)
        count.female <- apply(geno[sex==2,], 2, sum)
        n.male   <- sum(sex==1)
        n.female <- sum(sex==2)
        maf <- (count.male + count.female)/(n.male + 2*n.female)
      } else
    {
      maf <- apply(geno/2, 2, mean)
     }

    if(is.null(weights)) {
      wt <- 1/sqrt(maf * (1-maf))
    } else {
      wt <- weights[ v > 0 ]
    }

    # estimate cor among markers
    r.mat <- cor(geno)
   
    ## compute f from top of p412.  results in f=1 in MD weights,
    ## left here in case other weights used
    
    f <- wt * sqrt(maf*(1-maf))
    
    fRmat <- (f %o% f) * r.mat

    var.z <- fRmat * c.factor
    
    e.Q <- sum(diag(var.z))
    var.Q <- 2*sum(diag(var.z %*% var.z))
    
    # score males according to male.dose
    if(x.chrom & (male.dose !=1) )
      {
        geno.score <- geno
       
        geno.score[sex==1,] <- geno[sex==1,]*male.dose


        # kernel stat info
        kmat <- geno.score %*% diag(wt^2) %*% t(geno.score)

        # Burden stat info
        burden.score.subject <- as.vector(geno.score %*% wt)
        
      } else
    {
     
      # kernel stat info
      kmat <- geno %*% diag(wt^2) %*% t(geno)

      # Burden stat info
      burden.score.subject <- as.vector(geno %*% wt)
    }

   ########## Quadratic kernel stat
    
    stat.kernel <- as.vector(resid %*% kmat %*% resid)
    
    scale <- 2*e.Q / var.Q
    df    <- 2*e.Q^2 / var.Q
    pval.kernel <- pchisq(stat.kernel*scale, df=df, lower.tail=FALSE)
 
    eig <-eigen(var.z, symmetric=T, only.values=T)
    evals <-eig$values[eig$values>1e-6*eig$values[1]] 
    pval.kernel.davies <- davies(stat.kernel , evals, acc=accuracy.davies)$Qq

    ########### Burden stat (2-sided)
     
    stat.num <- (resid %*% burden.score.subject)^2
    factor.sum <- sum(fRmat)    
    stat.denom <- factor.sum * c.factor
    
    burden.stat <- stat.num / stat.denom
    burden.pval <- pchisq(burden.stat, 1,lower.tail=FALSE)
    lst <- list(stat.kernel = stat.kernel, df.kernel=df,
                     pval.kernel = pval.kernel,
		     pval.kernel.davies=pval.kernel.davies,
                     stat.burden = burden.stat,
                     pval.burden = burden.pval)

    return(lst)
    
}
