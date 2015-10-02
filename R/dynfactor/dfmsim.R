dfmsim <- function(L,Sigma,T,N,K){
#
#	dfmsim.R
#	Simulates data from a dynamic factor model
#
#	Original code was written by Dimitris Korobilis and adapted by J.C. Martinez-Ovando (2015)
#
#	Author: 	Juan~Carlos Martínez-Ovando
#	Email:		juan.martinez@banxico.org.mx
#				jc.martinez.ovando@gmail.com
#

if(nargs() == 0){
    T <- 200
    N <- 9
    K <- 3
    lag_f <- 1

    L <- t(matrix(c(1,0,0,
					0, 1, 0,
					0, 0, 1,
					0.99, 0.65, 0.34,
					0.99, 0, 0,
					0.78, 0, 0.54,
					0.35, 0.85, 0.78,
					0, 0.33, 0.90,
					0.78, 0, 0.90), 3, 9))

    Sigma <- diag(c(0.2, 0.19, 0.36, 0.2, 0.2, 0.19, 0.19, 0.36, 0.36))
    
    PHI <- t(matrix(c(0.5, 0.0, 0.0,
					  0.0, 0.5, 0.0,
					  0.0, 0.0, 0.5), 3, 3))
    
    PSI <- t(matrix(c(1.0, 0.5, 0.5, 
					  0.0, 1.0, 0.5,
					  0.0, 0.0, 1.0), 3, 3))
    
    Q <- solve(PSI %*% t(PSI))
	}

f  <- rbind(matrix(runif(lag_f*K, min = 0, max = 1), lag_f, K),  matrix(0, T, K))
#	Now generate f from VAR (L,PHI,PSI)
nn <- 1
for(nn in lag_f:(T+lag_f)){
    u <- t(chol(Q)) %*% matrix(runif(K),K,1)
    flag <- mlag(f,lag_f)
    f[nn, ] <- flag[nn, ] %*% PHI + t(u)
	}
f <- f[ (lag_f+1):(T+lag_f), ]
X <- f %*% t(L) + matrix(runif(T*N, min = 0, max = 1), T, N) %*% t(chol(Sigma))

return(X)
}
#
# -- End of "dfmsim.R"