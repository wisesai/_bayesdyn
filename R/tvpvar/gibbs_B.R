gibbs_B <- function(y,Z,Ht,Qdraw,K,M,t,B_0_prmean,B_0_prvar,Q_prmean,Q_prvar){
#
#	gibbs_B.R
#
#	Author: 	Juan~Carlos Martínez-Ovando
#	Email:		juan.martinez@banxico.org.mx
#				jc.martinez.ovando@gmail.com
#

#	Btdrawc is a draw of the mean VAR coefficients, B(t)
#
y1 <- y
Z1 <- Z
Ht1 <- Ht
Qt1 <- Qdraw
m1 <- K
p1 <- M
t1 <- t
B01 <- B_0_prmean
V01 <- B_0_prvar

#	Carte and Kohn
gibbs_ssa_out <- gibbs_ssa(y1, Z1, Ht1, Qt1, m1, p1, t1, B01, V01)
Btdrawc <- gibbs_ssa_out[[1]]
log_lik <- gibbs_ssa_out[[2]] 
	
#	Accept draw
Btdraw <- Btdrawc
#	or use the code below to check for stationarity
#     #	Now check for the polynomial roots to see if explosive
#     ctemp1 <- matrix(0, M, M*p_lags)
#     counter <- NULL
#     restviol <- 0
#	  i <- 1
#     for(i in 1:t){
#         BBtempor <- Btdrawc[ ,i]
#         BBtempor <- t(reshaper(BBtempor, 'col', M*p_lags, M))
#         ctemp1 <- rbind( BBtempor, diag(M*(p_lags-1)), matrix(0, M*(p_lags-1), M))
#         if(max(abs(eig(ctemp1)))>0.9999){
#             restviol <- 1
#             counter <- cbind(counter, restviol)
#         	}
#     	}
#     #	if they have been rejected keep old draw, otherwise accept new draw 
#     if(sum(counter)==0){
#         Btdraw <- Btdrawc;
#         print('I found a keeper!')
#     	}

#	Sampling Q -- the covariance of B(t) (from iWishart).
#	Take the SSE in the state equation of B(t).
Btdraw_aux1 <- as.matrix(t(Btdraw[ , 2:t]))
Btdraw_aux2 <- as.matrix(t(Btdraw[ , 1:(t-1)]))
Qdraw <- gibbs_Q(Btdraw_aux1,Btdraw_aux2,K,Q_prmean,Q_prvar)

#	Output
output <- list(Btdraw,Qdraw,log_lik)
return(output)
}
#
#	-- End of "gibbs_B.R" --