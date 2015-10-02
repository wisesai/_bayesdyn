wishartrnd <- function(h,n){
#
#	wishartrnd.R
#	This function generates random samples from a Wishart distribution
#	with 'h' degrees of freedom.
#	Purpose:	Draws an m x m matrix from a wishart distribution
#   	        with scale matrix h and degrees of freedom nu = n.
#   	        This procedure uses Bartlett's decomposition.
#	Inputs:		h     -- m x m scale matrix.
#   	        n     -- scalar degrees of freedom.
#	Outputs:	s     -- m x m matrix draw from the wishart
#
#	Author:
# 	Juan~Carlos Martínez-Ovando
#	Email:		juan.martinez@banxico.org.mx
#				jc.martinez.ovando@gmail.com
#   	                 distribution.
#
#	Note:		Parameterized so that mean is n*h
#

library('mvtnorm')

p.order <- dim(h)[1]
A <- t(chol(h)) %*% t(rmvnorm(n, mean = matrix(0,p.order), sigma = diag(p.order), method="chol"))
A <- A %*% t(A)

#	Output
return(A)
}
#
#	--	End of "wishartrnd.R" --