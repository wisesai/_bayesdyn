matrixcol <- function(M){
#
#	matrixcol.R
#	Convert the matrix M (dim 'n x p') into a comun vector (dim 'np x 1')
#	Conversion is made column by column in M
#
#	Author: 	Juan~Carlos Martínez-Ovando
#	Email:		juan.martinez@banxico.org.mx
#				jc.martinez.ovando@gmail.com
#

n <- dim(M)[1]
p <- dim(M)[2]

vec <- matrix(NaN,n*p,1)
k <- 1
cont <- 0
for(k in 1:p){
  vec[((n*cont)+1):((n*cont)+n)] <- M[, k]
  cont <- cont + 1
  }

#	Output
return(vec)
}
#
#	-- END of "matrixcol.R" --