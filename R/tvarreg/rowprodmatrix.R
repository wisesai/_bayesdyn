rowprodmatrix <- function(w,alpha_sim){
#
#	This function computes the row-by-row product between two matrices of the same dimension
#
#	Revision:		May 17, 2015
#

#	Dimensions
T <- dim(w)[1]
k <- dim(w)[2]

rowprod <- NaN * matrix(1, T, 1)
tt <- 1
for(tt in 1:T){
	rowprod[tt,1] <- t(w[tt, ]) %*% alpha_sim[tt, ]
  }

#	Output
return(rowprod)
}
#
#	-- End of "rowprodmatrix.R" --