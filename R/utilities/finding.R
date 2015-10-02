finding <- function(A,j){
#
#	finding.R
#	Finds the indexes of an array that fulfill an equality condition with j
#
#	Author: 	Juan~Carlos Martínez-Ovando
#	Email:		juan.martinez@banxico.org.mx
#				jc.martinez.ovando@gmail.com
#

if(any(A==j)==FALSE){
	indexes <- c(0)
}else{
	indexes <- which(A == j)
}

#	Output
return(indexes)
}
#
# -- End of "finding.R"