gibbs_ssb <- function(y,Z,Ht,Qt,m,p,t,B0,V0){
#
#	gibbs_ssb.R
#	This function implements the Gibbs sampler for state space models based on
#	the algorithm proposed by Carter and Kohn (1994)
#
#	Original code was written by Piotr Eliasz (2003) and modified by D. Korobilis (2007)
#	This code was written by J.C. Martinez-Ovando (2015)
#
#	Author: 	Juan~Carlos Martínez-Ovando
#	Email:		juan.martinez@banxico.org.mx
#				jc.martinez.ovando@gmail.com
#

#	Kalman Filter
bp <- B0
Vp <- V0
bt <- matrix(0, t, m)
Vt <- matrix(0, m^2, t)
log_lik <- 0
i <- 1
for(i in 1:t){
    R <- as.matrix(Ht[((i-1)*p+1):(i*p), ])
    H <- Z[((i-1)*p+1):(i*p), ]
    cfe <- as.matrix(y[ ,i]) - H %*% bp			#	conditional forecast error
    f <- t(H) %*% Vp %*% H + R					#	variance of the conditional forecast error
    inv_f <- solve(f)
    log_lik <- log_lik + log(det(f)) + t(cfe) %*% inv_f %*% cfe
    btt <- bp + Vp %*% H %*% inv_f %*% cfe
    Vtt <- Vp - Vp %*% H %*% inv_f %*% t(H) %*% t(Vp)
    if(i < t){
        bp <- btt
        Vp <- Vtt + Qt
		}
    bt[i, ] = t(btt)
    Vt[ ,i] <- reshaper(Vtt,'col')
	}

#	draw Sdraw(T|T) ~ N(S(T|T),P(T|T))
bdraw <- matrix(0, t, m)
bdraw[t, ] <- rmvnorm(1,btt,Vtt)

#	Backward recurssions
i <- 1
for(i in 1:(t-1)){
    bf <- t(bdraw[t-i+1, ])
    btt <- t(bt[t-i, ])
	Vtt <- reshapes(as.matrix(Vt[ ,t-i]), 'col', m, m)
    f <- Vtt + Qt
    inv_f <- solve(f)
    cfe <- bf - btt
    bmean <- t(btt) + Vtt %*% inv_f %*% t(cfe)
    bvar <- Vtt - Vtt %*% inv_f %*% Vtt
    bdraw[t-i, ] <- rmvnorm(1,bmean,bvar)
	}
bdraw <- t(bdraw)

#	output
output <- list(bdraw,log_lik)
return(output)
}
#
#	-- End of "gibbs_ssb.R" --