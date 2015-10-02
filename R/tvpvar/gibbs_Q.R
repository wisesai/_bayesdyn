gibbs_Q <- function(Btdraw_aux1,Btdraw_aux2,K,Q_prmean,Q_prvar){
#
#	gibbs_Q.R
#
#	Original code was written by Piotr Eliasz (2003) and modified by D. Korobilis (2007)
#	This code was written by J.C. Martinez-Ovando (2015)
#
#	Author: 	Juan~Carlos Martínez-Ovando
#	Email:		juan.martinez@banxico.org.mx
#				jc.martinez.ovando@gmail.com
#

Btemp <-  Btdraw_aux1 - Btdraw_aux2
sse_2 <- matrix(0, K, K)
i <- 1 
for(i in 1:(t-1)){
    sse_2 <- sse_2 + as.matrix(Btemp[i,]) %*% t(as.matrix(Btemp[i,]))
  }

#	Subsequently draw Q, the covariance matrix of B(t)
Qinv <- solve(sse_2 + Q_prmean)
Qinvdraw <- wishartrnd(Qinv, t+Q_prvar)
Qdraw <- solve(Qinvdraw)				#	this is a draw from Q

#	Output
return(Qdraw)
}
#
#	-- End of "gibbs_Q.R" --