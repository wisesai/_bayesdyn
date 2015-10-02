dynfactor <- function(X, demean=1, K, lags=1, N_mcmc, N_burn, N_thin=1, N_print){
#
#	dynfactor.R
#	This implements the state-space estimation of the dynamic factor model,
#	their loadings and autoregressive parameters. 
#
#  --------------------------------------------------------------------------
#   Estimate dynamic factor model:
#                  X_t = Lam x F_t + e_t
#                  F_t = B(L) x F_t-1 + u_t
#   The previous (Eliasz, 2005) version was:
#           [X_t ; Y_t] = [Lf Ly ; 0 I] x [F_t ; Y_t] + e_t
#           [F_t ; Y_t] = [B(L) ; 0 I] x [F_t-1 ; Y_t-1] + u_t
#  
#  --------------------------------------------------------------------------
#
#	Input:
#		X			-	
#		demean		-	
#		K			-	
#		lags		-	
#		N_mcmc		-	
#		N_burn			-	
#		N_thin		-	
#		N_print	-
#
#	Output list:
#		Ldraw		-	
#		Bdraw		-	
#		Qdraw		-	
#		Fdraw		-	
#		Rdraw		-	
#
#	Original code was written by Piotr Eliasz (2003) and modified by D. Korobilis (2007)
#	This code was written by J.C. Martinez-Ovando (2015)
#
#	Author: 	Juan~Carlos Martínez-Ovando
#	Email:		juan.martinez@banxico.org.mx
#				jc.martinez.ovando@gmail.com
#
#	Reference:	...
#

if(demean == 1){
	T <- dim(X)[1]
	N <- dim(X)[2]
	#	Demean xraw
	X <- X - matrix(1, T, N) %*% diag(colMeans(X))
  }

if(lags != 1){
	stop("dynfactor: is curently build for lags==1 only...")
  }
  
#	Store draws in:
Ldraw <- array(0, c(N_mcmc-N_burn, N, K) )
Bdraw <- array(0, c(N_mcmc-N_burn, K, K, lags))
Qdraw <- array(0, c(N_mcmc-N_burn, K, K) )
Fdraw <- array(0, c(N_mcmc-N_burn, T, K) )
Rdraw <- array(0, c(N_mcmc-N_burn, K^2) )

#	********************************************************
#	STANDARDIZE for PC only
X_st <- X / matrix(1, T, N) %*% diag(apply(X, 2, sd))
  
#	First step - extract Principal Components from X
XX <- X_st
k <- K
X_st_extract <- pcext(X_st,K)
F0 <- X_st_extract[[1]]
Lf <- X_st_extract[[2]] 

#	Transform factors and loadings for LE normalization
Lf_qr <- qr(t(Lf))
rl <- Lf_qr$qr
ql <- matrix(Lf_qr$qraux, sqrt(length(Lf_qr$qraux)), sqrt(length(Lf_qr$qraux)))

Lf <- rl								#	 do not transpose yet, is upper triangular
F0 <- F0 %*% ql

#	Need identity in the first K columns of Lf, call them A for now
A <- Lf[ ,1:K]
Lf <- t(cbind(diag(K), solve(A) %*% Lf[ ,(K+1):N]))
F0 <- F0 %*% A

#	Obtain R
e <- X_st - F0%*%t(Lf)
R <- t(e)%*%e / T
R <- diag(diag(R))
L <- Lf

#	Run a VAR in F, obtain initial B and Q
varfreq_out <- varfreq(F0,lags,NULL)
B <- varfreq_out[[1]]
Bc <- varfreq_out[[2]] 
v <- varfreq_out[[3]]
Q <- varfreq_out[[4]]
invFF <- varfreq_out[[5]]



#	Put it all in state-space representation, write obs equ as XY=FY*L+e
XY <- X									#	Tx(N+M)
FY <- F0

#	adjust for lags in state equation, Q is KxK
if(lags > 1){
	Q <- cbind(Q, matrix(0, K, K*(lags-1)))
	Q <- rbind(Q, matrix(0, K*(lags-1),K*lags))
	B <- rbind(B, diag(K*(lags-1)))
	B <- cbind(B, matrix(0, 2*K*(lags-1), K))
  }else{
	Q <- Q
	B <- B
  }

#	start with
Sp <- matrix(0, K*lags, 1)
Pp <- diag(K*lags)  

#	Proper Priors:-----------------------------------------------------------
#	on VAR -- Normal-Wishart, after Kadiyala, Karlsson, 1997
#	on Q -- si
#	on B -- increasing tightness
#	on observable equation:
#	N(0,I)-iG(3,0.001)

#	prior distributions for VAR part, need B and Q
vo <- K+2
s0 <- 3
alpha <- 0.001
L_var_prior <- diag(K)
Qi <- zeros(K,1)

#	singles out latent factors
indexnM <- matrix(1, K, lags)
indexnM <- which(indexnM==1)
#	***************End of Preliminaries & PriorSpecification******************

#	==============================================================================================
#	====================================== START SAMPLING ========================================
#	==============================================================================================
tic <- proc.time()																#	This is just a timer
print('Number of iterations')
mod <- function(x,m){t1<-floor(x/m); return(x-t1*m)}

irep <- 1
for(irep in 1:(N_mcmc + N_burn)){
	#	Print iterations
    if( (mod(irep,N_print) == 0) & (irep < N_burn) ){
	    print(paste("Burn-in: ", irep, " of ", N_burn, sep= ""))
	  }else if( (mod(irep,N_print) == 0) & (irep > N_burn) ){
	    print(paste("MCMC: ", irep-N_burn, " of ", N_mcmc, sep= ""))
	  }

    #	STEP 1. =========|DRAW FACTORS
    #	generate Gibbs draws of the factors
    H <- L
    F <- B
    t <- dim(XY)[1]
	n <- dim(XY)[2]
    kml <- dim(Sp)[1]
    km <- dim(L)[2]
    S <- matrix(0, t, kml)
    P <- matrix(0, kml^2, t)
    Sdraw <- matrix(0, t, kml)
	i <- 1
    for(i in 1:t){	# -- Start of "for 1"
        y <- as.matrix(XY[i,])
        nu <- y - H %*% Sp[1:km]				#	conditional forecast error
        f <- H%*%Pp[1:km,1:km]%*%t(H) + R		#	variance of the conditional forecast error
        finv <- solve(f)
        
        Stt <- Sp + Pp[ ,1:km]%*%t(H)%*%finv%*%nu
        Ptt <- Pp - Pp[ ,1:km]%*%t(H)%*%finv%*%H%*%Pp[1:km, ]
        
        if(i < t){
            Sp <- F%*%Stt
            Pp <- F%*%Ptt%*%t(F) + Q
          }
        
        S[i, ] <- t(Stt)
        P[ ,i] <- reshaper(Ptt, 'col')
	  } # -- End of "for 1" 
    
    #	draw Sdraw(T|T) ~ N(S(T|T),P(T|T))
    Sdraw[t, ] <- S[t, ]
    Sdraw[t, indexnM] <- mvrnorm(1, Sdraw[t,indexnM],Ptt[indexnM,indexnM])

    #	iterate 'down', drawing at each step, use modification for singular Q
    Qstar <- Q[1:km, 1:km]
    Fstar <- F[1:km, ]	

    i <- (t-1)
	for(i in (t-1):1){
        Sf <- as.matrix(Sdraw[i+1, 1:km])
        Stt <- as.matrix(S[i, ])
        Ptt <- reshapes(as.matrix(P[ ,i]),'col', kml, kml)
        f <- Fstar%*%Ptt%*%t(Fstar) + Qstar
        finv <- solve(f)
        nu <- Sf - Fstar%*%Stt
        
        Smean <- Stt + Ptt%*%t(Fstar)%*%finv%*%nu
        Svar <- Ptt - Ptt%*%t(Fstar)%*%finv%*%Fstar%*%Ptt
       
        Sdraw[i, ] <- t(Smean)
        Sdraw[i, indexnM] <- mvrnorm(1, Sdraw[i,indexnM], Svar[indexnM, indexnM])
      }
    FY <- Sdraw[ ,1:km]	

	#	Demean
    FY <- FY - matrix(1,T,dim(FY)[2]) %*% diag(colMeans(FY))

    #	STEP 2. ========|DRAW COEFFICIENTS
    #	-----------------------2.1. STATE EQUATION---------------------------
    #	first univ AR for scale in priors
    i <- 1
	Qi <- 0*c(1:km)
	for(i in 1:km){
		varfreq_out2 <- varfreq(as.matrix(FY[ ,i]),lags,NULL)
		Bi <- varfreq_out2[[1]]
		Bci <- varfreq_out2[[2]]
		vi <- varfreq_out2[[3]]
		Qi[i] <- varfreq_out2[[4]]
		invFYFYi <- varfreq_out2[[5]]
      }
    Q_prior <- diag(Qi)
    B_var_prior <- diag(kronecker(1/Qi,1/(1:lags)))
    estvar_out2a <- varfreq(FY,lags,NULL)
    Bd <- estvar_out2a[[1]]
	Bdc <- estvar_out2a[[2]]
	v <- estvar_out2a[[3]]
	Qd <- estvar_out2a[[4]]
	invFYFY <- estvar_out2a[[5]]
	B_hat <- t(Bd)
    Z <- array(0, c(T, km, lags))
	i <- 1
    for(i in 1:lags){
        Z[(lags+1):T, , i] <- FY[(lags+1-i):(T-i), ]
      }
    Z <- Z[ , , 1]
    Z <- Z[(lags+1):T, ]
	iB_var_prior <- solve(B_var_prior)
    B_var_post <- solve(iB_var_prior + t(Z)%*%Z)
    B_post <- B_var_post %*% (t(Z)%*%Z) %*% B_hat
    Q_post <- t(B_hat) %*% t(Z) %*% Z %*% B_hat + Q_prior + (T-lags)*Qd - t(B_post) %*% (iB_var_prior + t(Z)%*%Z) %*% B_post
    
    #	Draw Q from inverse Wishart
    iQd <- matrix(rnorm((T+vo)*km, mean = 0, sd = 1), T+vo,km) %*% chol(solve(Q_post))
    Qd <- solve(t(iQd)%*%iQd)
    Q[1:km,1:km] <- Qd

    #	Draw B conditional on Q
    vecB_post <- reshaper(B_post,'col')
    vecBd <- vecB_post + t(chol(kronecker(Qd,B_var_post))) %*% rnorm(km*km*lags, mean = 0, sd = 1)
    Bd <- t(reshapes(vecBd,'col',km*lags,km))
    B[1:km, ] <- Bd
    
    #	Truncate to ensure stationarity
    while( max(abs(eig(B))) > 0.999){
        vecBd <- vecB_post + t(chol(kronecker(Qd, B_var_post))) %*% rnorm(km*km*lags, mean = 0, sd = 1) 
		Bd <- t(reshapes(vecBd,'col',km*lags,km))
        B[1:km, ] <- Bd
      }
	
    #	----------------------2.2. OBSERVATION EQUATION----------------------
    #	OLS quantities
    L_OLS <- solve(t(FY)%*%FY)%*%(t(FY)%*%X[ ,(K+1):N])
    R_OLS <- t(X - FY%*%t(L)) %*% (X - FY%*%t(L)) * (T-N)^-1
    
    
    L <- t(cbind(diag(K), L_OLS))

	n <- 1
	for(n in 1:N){
        ed <- X[ ,n] - FY%*%as.matrix(L[n,])

        #	draw R(n,n)
        R_bar <- s0 + t(ed)%*%ed + t(as.matrix(L[n, ]))%*%solve(L_var_prior+solve(t(FY)%*%FY))%*%as.matrix(L[n, ])        
        Rd <- rchisq(1, T+alpha, ncp = 0)
        Rd <- R_bar/Rd
        R[n, n] <- Rd
      }

    #	Save draws
    if(irep > N_burn){
        Ldraw[(irep-N_burn), , ] <- L[1:N, 1:km]
        Bdraw[(irep-N_burn), , , ] <- Bd					#	reshaper[Bd,'col', km, km, lags]	#	in Matlab
        Qdraw[(irep-N_burn), , ] <- Qd
        Fdraw[(irep-N_burn), , ] <- FY
        Rdraw[(irep-N_burn), ] <- as.matrix(diag(R))
      }
	  
  }# -- End of "for 'irep'-- GIBBS iterations... "
#	====================== End Sampling Posteriors ===========================

print('Time-processing elapsed:')
proc.time() - tic

#	Output 
output <- list(Ldraw,Bdraw,Qdraw,Fdraw,Rdraw)
return(output)
}
#
# -- End of "dynfactor.R" --