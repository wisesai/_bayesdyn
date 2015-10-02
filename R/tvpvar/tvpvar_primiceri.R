tvpvar_primiceri <- function(Y,tau,M,p){
#
#	tvpvar_primiceri.R
#	Time Series prior under Primiceri's approach
#
#	Original code was written by D. Korobilis and adpated by J.C. Martinez-Ovando (2015)
#
#	Author: 	Juan~Carlos Martínez-Ovando
#	Email:		juan.martinez@banxico.org.mx
#				jc.martinez.ovando@gmail.com
#

yt <- t(Y[(p+1):(tau+p), ])

#	m is the number of elements in the state vector
m <- M + p*(M^2)
Zt <- NULL

i <- (p+1)
for(i in (p+1):(tau+p)){
    ztemM <- diag(M)
    j <- 1
	for(j in 1:p){
        xlag <- Y[i-j, 1:M]
        xtemM <- matrix(0, M, M*M)
		jj <- 1
        for(jj in 1:M){
            xtemM[jj, ((jj-1)*M+1):(jj*M)] <- xlag;
			} # -- End of "for 3"
        ztemM <- cbind(ztemM, xtemM)
		} # -- End of "for 2"
    Zt <- rbind(Zt,ztemM)
	} #	-- End of "for 1"

vbar <- matrix(0, m, m)
xhy <- matrix(0, m, 1)
i <- 1
for(i in 1:tau){
    zhat1 <- Zt[((i-1)*M+1):(i*M), ]
    vbar <- vbar + t(zhat1) %*% zhat1
    xhy <- xhy + t(zhat1) %*% yt[ ,i]
	}

vbar <- solve(vbar)
aols <- vbar %*% xhy

sse2 <- matrix(0, M, M)
i <- 1
for(i in 1:tau){
    zhat1 <- Zt[((i-1)*M+1):(i*M), ]
    sse2 <- sse2 + (yt[ ,i] - zhat1 %*% aols) %*% t(yt[ ,i] - zhat1 %*% aols)
	}
hbar <- sse2 / tau

vbar <- matrix(0, m, m)
i <- 1
for(i in 1:tau){
    zhat1 <- Zt[((i-1)*M+1):(i*M), ]
    vbar <- vbar + t(zhat1) %*% solve(hbar) %*% zhat1
	}

vbar <- solve(vbar)
achol <- chol(hbar)
ssig <- matrix(0, M, M)
i <- 1
for(i in 1:M){
    ssig[i, i] <- achol[i, i] 
	j <- 1
    for(j in 1:M){
        achol[j, i]  <- achol[j, i] / ssig[i,i]
		}# --  End of "for 2"
	}# --  End of "for 1"

achol <- solve(achol)
numa <- M*(M-1)/2
a0 <- matrix(0, numa, 1)
ic <- 1
i <- 2
for(i in 2:M){
	j <- 1
    for(j in 1:(i-1)){
        a0[ic, 1] <- achol[i, j]
        ic <- ic + 1
		} # -- End of "for 1"
	} # -- End of "for 1" 
ssig1 <- matrix(0, M, 1)
i <- 1
for(i in 1:M){
    ssig1[i, 1] <- log(ssig[i, i]^2)
	}

hbar1 <- solve(tau*hbar)
hdraw <- matrix(0, M, M)
a02mo <- matrix(0, numa, numa)
a0mean <- matrix(0, numa, 1)
ireM <- 1
for(ireM in 1:4000){
	h <- hbar1
	n <- tau
    hdraw <- wishartrnd(h,n)
    hdraw <- solve(hdraw)
    achol <- t(chol(hdraw))
    ssig <- matrix(0, M, M)
	i <- 1
    for(i in 1:M){
        ssig[i, i] <- achol[i,i]
		j <- 1
        for(j in 1:M){
            achol[j, i] <- achol[j, i] / ssig[i, i]
			} # -- End of "for 2" 
		}# -- End of "for 1"
    achol <- solve(achol)
    a0draw <- matrix(0, numa, 1)
    ic <- 1
	i <- 2
    for(i in 2:M){
		j <- 1
        for(j in 1:(i-1)){
            a0draw[ic, 1] <- achol[i, j]
            ic <- ic + 1;
			} # End of "for 2"
		} # -- End of "for 1"
    a02mo <- a02mo + a0draw %*% t(a0draw)
    a0mean <- a0mean + a0draw
	} # -- End of "for 1"

a02mo <- a02mo /4000
a0mean <- a0mean /4000
a02mo <- a02mo - a0mean %*% t(a0mean)

#	Output	
output <- list(aols,vbar,a0,ssig1,a02mo)
return(output)
}
#
#	-- End of "tvpvar_primiceri.R" --