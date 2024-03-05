
## Title: format.relation.R
## Purpose: worker function to reformat the relationship matrix to evaluate
## 			the kinship matrix when special relationships between pairs of 
##			individuals exists
## Author: Alessia Visconti
## Created: 20190411
## Updated: 20190411

format_relation <- function(relation)

##Arguments:
##  
##  Required:
##  relation:  data.frame/matrix with 4 columns (id, person1, person2, code) specifying 
##	special relationships between pairs of individuals and used by the kinship 
##	function. Codes are : 
##		1=Monozygotic twin, 
##		2=Dizygotic twin, 
##		3=Twin of unknown zygosity, and
##		4=Spouse. 

{
	relation <- as.data.frame(relation)
	names(relation) <- c("ped", "id1", "id2", "code")
	relation$id1 <- paste(relation$ped, relation$id1, sep="-")
	relation$id2 <- paste(relation$ped, relation$id2, sep="-")
	relation$ped <- NULL
	relation
}
