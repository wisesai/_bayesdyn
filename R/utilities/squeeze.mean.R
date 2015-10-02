squeeze.mean <- function(Array){
#
#	squeeze.sd.R
#	Computes the mean of a three-dimensional array (dim JxKxL) with respect to the first dimension.
#	Output is a matrix array of dimension KxL.
#
#	Author: 	Juan~Carlos Martínez-Ovando
#	Email:		juan.martinez@banxico.org.mx
#				jc.martinez.ovando@gmail.com
#

if(length(dim(Array)) != 3){stop("squeeze.sd::: Error in the specification of 'Array'.")}
J <- dim(Array)[1]
K <- dim(Array)[2]
L <- dim(Array)[3]

Array.Saqueeze <- matrix(NaN,K,L)
k <- 1
for(k in 1:K){
  l <- 1
  for(l in 1:L){
    Array.Saqueeze[k,l] <- sd(Array[,k,l])
	}
  }

#	Output
return(Array.Saqueeze)
}
#
# -- End of "squeeze.sd.R"