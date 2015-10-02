mlag2 <- function(Y,p){
#
#	mlag2.R
#	PURPOSE: generates a matrix of n lags from a matrix (or vector)
#           containing a set of vectors (For use in var routines)
#
#	Author:
# 	Juan~Carlos Martínez-Ovando
#	Email:		juan.martinez@banxico.org.mx
#				jc.martinez.ovando@gmail.com
#

Y <- as.matrix(Y)
Traw <- dim(Y)[1]
N <- dim(Y)[2]

Np <- N*p
Ylag <- matrix(0, Traw, Np)

ii <- 1
for(ii in 1:p){
    ylag_cols <- (N*(ii-1)+1):(N*ii)
	ylag_rows <- (p+1):Traw
	y_rows <- (p+1-ii):(Traw-ii)
	y_cols <- 1:N
	Ylag[ylag_rows, ylag_cols] <- Y[y_rows, y_cols]
	}

#	Output
return(Ylag)
}
#
# -- End of "mlag2.R" --