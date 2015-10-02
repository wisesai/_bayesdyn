indexownlags <- function(i,constant,M,p,K){
#
#	indexownlags.R
#	Index in each equation which are the own lags.
#
#	Author: 	Juan~Carlos Martínez-Ovando
#	Email:		juan.martinez@banxico.org.mx
#				jc.martinez.ovando@gmail.com
#

index <- matrix(NaN,1,p)

reference <- round(length(M:K)/2) + 1
index.aux <- c(M:K)

if( (constant+1) == 2 ){
  index[1,1] <- constant + i
  index[1,2] <- index.aux[reference + (i-1)]
  }else{
  index[1,1] <- constant + i
  index[1,2] <- index.aux[reference + (i-1)]}
 
#	Output 
return(index)	
}
#
# -- End of "indexownlags.R"