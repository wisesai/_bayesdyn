corrvc <- function(vc){
#
#	corrvc.R
#	Computes a correlation matrix from a variance-covariance matrix.
#              
#	Format:     cx <- corrvc(vc);
# 
#	Input:      vc    KxK variance-covariance matrix (of data or parameters).
# 
#	Output:     cx    KxK correlation matrix.
#
#	Author: 	Juan~Carlos Martínez-Ovando
#	Email:		juan.martinez@banxico.org.mx
#				jc.martinez.ovando@gmail.com
#

#	Check for complex input */
if(is.double(vc) == FALSE){       
    print('ERROR: Not implemented for complex arguments...')
}
std <- sqrt(diag(vc))
   
output <- vc / (std%*%t(std));

#	Output
return(output)
}
#
#	--	End of "corrvc.R"