colStdev <- function(x, na.rm = FALSE, dims = 1L){
#
#	colStdev.R
#	Computes the standard deviation of matrix 'x' by columns.
#
#	Author: 	Juan~Carlos Martínez-Ovando
#	Email:		juan.martinez@banxico.org.mx
#				jc.martinez.ovando@gmail.com
#

    if (is.data.frame(x)) 
        x <- as.matrix(x)
    if (!is.array(x) || length(dn <- dim(x)) < 2L) 
        stop("'x' must be an array of at least two dimensions")
    if (dims < 1L || dims > length(dn) - 1L) 
        stop("invalid 'dims'")
    n <- prod(dn[1L:dims])
    dn <- dn[-(1L:dims)]
    z <- matrix(NaN, 1, dn)
	if(is.complex(x)){ 
		k <- 1
		for(k in 1:dn){
			z[1,k] <- sd(Re(x[,k])) + (0+1i) * sd(Im(x[,k]))
		  }  
      }else{
		k <- 1
		for(k in 1:dn){
			z[1,k] <- sd(Re(x[,k]))
		  }  
	  }
#    if (length(dn) > 1L){
#       dim(z) <- dn
#        dimnames(z) <- dimnames(x)[-(1L:dims)]
#    }else{
#		names(z) <- dimnames(x)[[dims + 1]]
#	}

#	Output
return(t(z))
}
#
# -- End of "colStdev.R" --