squeeze.quantile <- function(Array,probs){
#
#	squeeze.quantile.R
#	Computes quantiles of a three-dimensional array (dim JxKxL).
#
#	Output:		Matrix array of dimension KxL
#
#	Author: 	Juan~Carlos Martínez-Ovando
#	Email:		juan.martinez@banxico.org.mx
#				jc.martinez.ovando@gmail.com
#

if(length(dim(Array)) != 3){print("squeeze.mean::: Error in the specification of 'Array'.")}
J <- dim(Array)[1]
K <- dim(Array)[2]
L <- dim(Array)[3]

Array.Saqueeze <- matrix(NaN,K,L)
k <- 1
for(k in 1:K){
  l <- 1
  for(l in 1:L){
    Array.Saqueeze[k,l] <- mean(Array[,k,l])
	}
  }

#	Output
return(Array.Saqueeze)
}
#
# -- End of "squeeze.quantile.R"