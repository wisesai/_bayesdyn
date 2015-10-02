pcext <- function(XX,k){
#
#	pcext.R
#	Extracts first k principal components from the matrix XX (dim 't x n'), 
#	loadings are normalized so that lam'lam/n=I, fac is t*k, lam is n*k
#
#	Original code was written by Dimitris Korobilis and adapted by J.C. Martinez-Ovando (2015)
#
#	Author: 	Juan~Carlos Martínez-Ovando
#	Email:		juan.martinez@banxico.org.mx
#				jc.martinez.ovando@gmail.com
#
t <- dim(XX)[1]
n <- dim(XX)[2]

xx <- t(XX) %*% XX

xx_eig <- eigen(xx)					#	eigvec correspond to eigval in descending order
eigval <- xx_eig[[1]]
eigvec <- xx_eig[[2]]

lam <- sqrt(n)*eigvec[ ,1:k]

fac <- XX %*% lam/n

output <- list(fac,lam)
return(output)
}
#
# -- End of "pcext.R"