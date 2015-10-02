mlag <- function(x, n, init){
#
#	mlag.R
#	PURPOSE: generates a matrix of n lags from a matrix (or vector)
#           containing a set of vectors (For use in var routines)
#	---------------------------------------------------
#	USAGE:     xlag = mlag(x,nlag)
#        or: xlag1 = mlag(x), which defaults to 1-lag
#	where: x = a matrix (or vector), nobs x nvar
#      nlag = # of contiguous lags for each vector in x
#      init = (optional) scalar value to feed initial missing values
#             (default = 0)
#	---------------------------------------------------
#	RETURNS:
#         xlag = a matrix of lags (nobs x nvar*nlag)
#         x1(t-1), x1(t-2), ... x1(t-nlag), x2(t-1), ... x2(t-nlag) ...
#  --------------------------------------------------
#	SEE ALSO: lag() which works more conventionally
#	---------------------------------------------------
#	Written in Matlab by:
#	James P. LeSage, Dept of Economics
#	University of Toledo
#	2801 W. Bancroft St,
#	Toledo, OH 43606
#	jpl@jpl.econ.utoledo.edu
#
#	Adapted to R by:
# 	Juan~Carlos Martínez-Ovando
#	Email:		juan.martinez@banxico.org.mx
#				jc.martinez.ovando@gmail.com
#

if(nargs() ==1){ 
	n <- 1							#  default value
	init <- 0
  }else if(nargs() == 2){
	init <- 0
  }

if(nargs() > 3){
  print('mlag: Wrong # of input arguments')
  }

nobs <- dim(x)[1]
nvar <- dim(x)[2]

xlag <- matrix(0, nobs, nvar*n)
icnt <- 0
i <- 1
for(i in 1:nvar){
	j <- 1
	for(j in 1:n){
		xlag[(j+1):nobs, icnt+j] <- x[1:(nobs-j), i]
		}
	icnt <- icnt + n
	}

#	Output
return(xlag)
}
#
# -- End of "mlag.R"