tvpvar <- function(ydata, demean=0, N_training_prior, p_lags, N_mcmc, N_burn, N_print, prior_type=1){
#
#	tvpvar.R
#	TVP-VAR Time varying structural VAR with stochastic volatility
#	------------------------------------------------------------------------------------
#	The model is:
#
#     Y(t) = B0(t) + B1(t)xY(t-1) + B2(t)xY(t-2) + e(t) 
# 
#	with e(t) ~ N(0,SIGMA(t)), and  L(t)' x SIGMA(t) x L(t) = D(t)*D(t),
#             _                                          _
#            |    1         0        0       ...       0  |
#            |  L21(t)      1        0       ...       0  |
#    L(t) =  |  L31(t)     L32(t)    1       ...       0  |
#            |   ...        ...     ...      ...      ... |
#            |_ LN1(t)      ...     ...    LN(N-1)(t)  1 _|
# 
# 
#	and D(t) = diag[exp(0.5 x h1(t)), .... ,exp(0.5 x hn(t))].
#
#	The state equations are
#
#            B(t) = B(t-1) + u(t),            u(t) ~ N(0,Q)
#            l(t) = l(t-1) + zeta(t),      zeta(t) ~ N(0,S)
#            h(t) = h(t-1) + eta(t),        eta(t) ~ N(0,W)
#
#	where 
#			B(t) 	=	[B0(t),B1(t),B2(t)]', 
#			l(t)	=	[L21(t),...,LN(N-1)(t)]'
#			h(t)	=	[h1(t),...,hn(t)]'.
#
#	Original code was written by D. Korobilis and adpated by J.C. Martinez-Ovando (2015)
#
#	Author: 	Juan~Carlos Martínez-Ovando
#	Email:		juan.martinez@banxico.org.mx
#				jc.martinez.ovando@gmail.com
#
#	************************************************************************************
#   NOTE: 
#      There are references to equations of Primiceri, "Time Varying Structural Vector
#      Autoregressions & Monetary Policy",(2005),Review of Economic Studies 72,821-852.
#	------------------------------------------------------------------------------------

#	---------------------------------------------------------------------------------
#	----------------------------	Preliminaries	---------------------------------
#	Demean and standardize data
if(demean == 1){
	t2 <- dim(ydata)[1]
	ydata_means <- t(as.matrix(colMeans(ydata)))
	ydata_stdev <- t(as.matrix(colStdev(ydata)))
	Y <- ( ydata- matrix(1,t2,1) %*% ydata_means ) / matrix(1,t2,1) %*% ydata_stdev
  }else{
	Y <- ydata
  }

#	Number of observations and dimension of X and Y
t <- dim(Y)[1]					# 	t is the time-series observations of Y
M <- dim(Y)[2]					#	M is the dimensionality of Y

numa <- M*(M-1)/2				# 	Number of lower triangular elements of A_t (other than 0's and 1's)

#	----------------------------	VAR equations	---------------------------------
#	Generate lagged Y matrix. This will be part of the X matrix
ylag  <- mlag2(Y,p_lags) 		# 	Y is [T x M]. ylag is [T x (Mp)]
#	Form RHS matrix 
ylag <- ylag[(p_lags+N_training_prior+1):t, ]	#	X_t = [1 y_t-1 y_t-2 ... y_t-k] for t=1:T

K <- M + p_lags*(M^2)			#	K is the number of elements in the state vector
# Create Z_t matrix.
Z <- matrix(0, (t-N_training_prior-p_lags)*M, K)

i <- 1
for(i in 1:(t-N_training_prior-p_lags)){
    ztemp <- diag(M)
	j <- 1
    for(j in 1:p_lags){        
        xtemp <- ylag[i, ((j-1)*M+1):(j*M)]
        xtemp <- t(kronecker(diag(M), xtemp))
        ztemp <- cbind(ztemp, xtemp)
		} # -- End "for 2"
    Z[((i-1)*M+1):(i*M), ] <- ztemp
	} # -- End "for 1"

#	Redefine TVP-VAR variables y
y <- t(Y[ (N_training_prior+p_lags+1):t, ])

#	Time series observations
t <- dim(y)[2]					#	t = (T - p_lags - N_training_prior)

#	--------------------------	Prior specification	-------------------------------
tvpvar_prior_out <- tvpvar_prior(prior_type,Y,N_training_prior,M,p_lags)
B_OLS <- tvpvar_prior_out[[1]]
VB_OLS <- tvpvar_prior_out[[2]]
A_OLS <- tvpvar_prior_out[[3]]
VA_OLS <- tvpvar_prior_out[[4]]
sigma_OLS <- tvpvar_prior_out[[5]]

# 	Set some hyperparameters
k_Q <- 0.01
k_S <- 0.1
k_W <- 0.01

#	Sizing some matrices as prior hyperparameters 
sizeW <- M					#	Size of matrix W
sizeS <- 1:M				#	Size of matrix S

#	Set Kalman filter initial conditions for the time-varying
#	parameters B(t), A(t) and log-Sigma(t). These are the mean VAR
#	coefficients, the lower-triangular VAR covariances and the diagonal
#	log-volatilities, respectively 

#	B_0 ~ N(B_OLS, 4Var(B_OLS))
B_0_prmean <- B_OLS
B_0_prvar <- 4*VB_OLS

#	A_0 ~ N(A_OLS, 4Var(A_OLS))
A_0_prmean <- A_OLS
A_0_prvar <- 4*VA_OLS

# log(Sigma_0) ~ N(log(sigma_OLS),I_n)
sigma_prmean <- sigma_OLS
sigma_prvar <- 4*diag(M)

#	Note that for IW distribution the _prmean/_prvar notation is preserved:
#	Q - is the covariance of B(t), 
#	S - is the covariance of A(t), and 
#	W - is the covariance of log-Sigma(t)

#	Q ~ IW(k2_Q*size(subsample)*Var(B_OLS),size(subsample))
Q_prmean <- ((k_Q)^2) * N_training_prior * VB_OLS
Q_prvar <- N_training_prior

#	W ~ IW(k2_W*(1+dimension(W))*I_n,(1+dimension(W)))
W_prmean <- ((k_W)^2) * (1 + sizeW) * diag(M)
W_prvar <- 1 + sizeW

# S ~ IW(k2_S*(1+dimension(S)*Var(A_OLS),(1+dimension(S)))
S_prmean <- list()		#	To be of size (M-1)x1	
S_prvar <- matrix(0, M-1, 1)
ind <- 1
ii <- 2
for(ii in 2:M){
    #	S is block diagonal as in Primiceri (2005)
    S_prmean[[ii-1]] <- ((k_S)^2)*(1 + sizeS[ii-1]) * VA_OLS[((ii-1)+(ii-3)*(ii-2)/2):ind,((ii-1)+(ii-3)*(ii-2)/2):ind]
    S_prvar[ii-1] <- 1 + sizeS[ii-1]
    ind <- ind + ii
	}

#	Parameters of the 7 component mixture approximation to a log(chi^2) density:
q_s <- c(0.00730, 0.10556, 0.00002, 0.04395, 0.34001, 0.24566, 0.25750)			#	probabilities
q_s <- as.matrix(q_s, length(q_s), 1)
m_s <- c(-10.12999, -3.97281, -8.56686, 2.77786, 0.61942, 1.79518, -1.08819)	#	means
m_s <- as.matrix(m_s, length(m_s), 1)
u2_s <- c(5.79596, 2.61369, 5.17950, 0.16735, 0.64009, 0.34023, 1.26261)		#	variances
u2_s <- as.matrix(u2_s, length(u2_s), 1)

#	--------------------------	Intialize matrices	-------------------------------
#	Specify covariance matrices for measurement and state equations
consQ <- 0.0001
consS <- 0.0001
consH <- 0.01
consW <- 0.0001
Ht <- kronecker(matrix(1, t, 1),consH*diag(M))				#	Initialize Htdraw, a draw from the VAR covariance matrix
Htchol <- kronecker(matrix(1, t, 1),sqrt(consH)*diag(M))	#	Cholesky of Htdraw defined above
Qdraw <- consQ*diag(K)										#	Initialize Qdraw, a draw from the covariance matrix Q
Sdraw <- consS*diag(numa)									#	Initialize Sdraw, a draw from the covariance matrix S
Sblockdraw <- list()										#	#cell(M-1,1)	#	...and then get the blocks of this matrix (see Primiceri)
ijc <- 1
jj <- 2
for(jj in 2:M){
    Sblockdraw[[jj-1]] <- Sdraw[((jj-1)+(jj-3)*(jj-2)/2):ijc,((jj-1)+(jj-3)*(jj-2)/2):ijc]
    ijc <- ijc + jj
	}
Wdraw <- consW*diag(M)										#	Initialize Wdraw, a draw from the covariance matrix W
Btdraw <- matrix(0, K, t)									#	Initialize Btdraw, a draw of the mean VAR coefficients, B(t)
Atdraw <- matrix(0, numa, t)								#	Initialize Atdraw, a draw of the non 0 or 1 elements of A(t)
Sigtdraw <- matrix(0, M, t)									#	Initialize Sigtdraw, a draw of the log-diagonal of SIGMA(t)
sigt <- kronecker(matrix(1, t, 1),0.01*diag(M))				#	Matrix of the exponent of Sigtdraws (SIGMA(t))
statedraw <- 5*matrix(1, t, M)								#	initialize the draw of the indicator variable 
															#	(of 7-component mixture of Normals approximation)
Zs <- kronecker(matrix(1, t, 1),diag(M))
prw <- matrix(0, length(q_s), 1)

#	Storage matrices for posteriors and stuff
Bt_postmean <- matrix(0, K, t)								#	regression coefficients B(t)
At_postmean <- matrix(0, numa, t)							#	lower triangular matrix A(t)
Sigt_postmean <- matrix(0, M, t)							#	diagonal std matrix SIGMA(t)
Qmean <- matrix(0, K, K)									#	covariance matrix Q of B(t)
Smean <- matrix(0, numa, numa)								#	covariance matrix S of A(t)
Wmean <- matrix(0, M, M)									#	covariance matrix W of SIGMA(t)

sigmean <- matrix(0, t, M)									#	mean of the diagonal of the VAR covariance matrix
cormean <- matrix(0, t, numa)								#	mean of the off-diagonal elements of the VAR cov matrix
sig2mo <- matrix(0, t, M)									#	squares of the diagonal of the VAR covariance matrix
cor2mo <- matrix(0, t, numa)								#	squares of the off-diagonal elements of the VAR cov matrix

#	-------------------------	End of preliminaries	-----------------------------
#	---------------------------------------------------------------------------------

#	---------------------------------------------------------------------------------
#	-----------------------------	MCMC sampling	---------------------------------
tic <- proc.time()											#	This is just a timer
print('Number of iterations')
mod <- function(x,m){t1<-floor(x/m); return(x-t1*m)}

irep <- 1
for(irep in 1:(N_mcmc + N_burn)){							#	Gibbs sampler iterations starts here
	#	Print iterations
    if( (mod(irep,N_print) == 0) & (irep < N_burn) ){
	    print(paste("Burn-in: ", irep, " of ", N_burn, sep= ""))
	  }else if( (mod(irep,N_print) == 0) & (irep > N_burn) ){
	    print(paste("MCMC: ", irep-N_burn, " of ", N_mcmc, sep= ""))
	  }
	  
    #	-----------------------------------------------------------------------------------------
    #		Step I:		Sampling B from p_lags(B|y,A,Sigma,V)
    #	-----------------------------------------------------------------------------------------
	#
	gibbs_B_out <- gibbs_B(y,Z,Ht,Qdraw,K,M,t,B_0_prmean,B_0_prvar,Q_prmean,Q_prvar)
	Btdraw <- gibbs_B_out[[1]]
	Qdraw <- gibbs_B_out[[2]]
	log_lik <- gibbs_B_out[[3]]

    #	-------------------------------------------------------------------------------------------
    #   	Step II: 	Sampling A(t) from p_lags(At|y,B,Sigma,V)
    #	-------------------------------------------------------------------------------------------
    #
	gibbs_A_out <- gibbs_A(M,t,y,Z,Btdraw,sigt,Sblockdraw,A_0_prmean,A_0_prvar,numa,S_prmean,S_prvar)
	Atdraw <- gibbs_A_out[[1]]
	Sinvblockdraw <- gibbs_A_out[[2]]
	yhat <- gibbs_A_out[[3]]
	
    #	------------------------------------------------------------------------------------------
    #   	STEP III:	Sampling diagonal VAR covariance matrix log-SIGMA(t)
    #	------------------------------------------------------------------------------------------
    #
	gibbs_logSigma_out <- gibbs_logSigma(M,t,Atdraw,yhat,statedraw,u2_s,m_s,Zs,Wdraw,sigma_prmean,sigma_prvar,q_s,W_prmean,W_prvar)
	Wdraw <- gibbs_logSigma_out[[1]]
	Ht <- gibbs_logSigma_out[[2]]
	Htsd <- gibbs_logSigma_out[[3]]
	
    #	----------------------------SAVE AFTER-BURN-IN DRAWS -----------------
    if(irep > N_burn){               
        #	Save only the means of parameters. Not memory efficient to store all draws (at least 
		#	for the time-varying parameters vectors, which are large). 
		#
		#	Note: If you want to store all draws, it is better to save them in a file at each iteration. 
        Bt_postmean <- Bt_postmean + Btdraw						#	regression coefficients B(t)
        At_postmean <- At_postmean + Atdraw						#	lower triangular matrix A(t)
        Sigt_postmean <- Sigt_postmean + Sigtdraw				#	diagonal std matrix SIGMA(t)
        Qmean <- Qmean + Qdraw									#	covariance matrix Q of B(t)
        ikc <- 1
		kk <- 2
        for(kk in 2:M){
            Sdraw[(((kk-1)+(kk-3)*(kk-2)/2)):ikc,(((kk-1)+(kk-3)*(kk-2)/2)):ikc] <- Sblockdraw[[kk-1]]
            ikc <- ikc + kk
			} # -- End of "for"
        Smean <- Smean + Sdraw									#	covariance matrix S of A(t)
        Wmean <- Wmean + Wdraw									# 	covariance matrix W of SIGMA(t)
        #	Get time-varying correlations and variances
        stemp6 <- matrix(0, 1, M)
        stemp5 <- matrix(0, t, M)
        stemp7 <- matrix(0, t, M)
		i <- 1
        for(i in 1:t){
            stemp8 <- corrvc(Ht[((i-1)*M+1):(i*M),])
            stemp7a <- NULL
            ic <- 1
			j <- 1
            for(j in 1:M){
                if(j > 1){
                    stemp7a <- rbind(stemp7a, as.matrix(stemp8[j,1:ic]))
                    ic <- ic + 1
				  }
                stemp6[1, j] <- sqrt(Ht[((i-1)*M+j),j])
				} # -- End of "for"
			stemp5[i, ] <- stemp6
            stemp7[i, ] <- stemp7a
			} # -- End of "for"
		stemp5 <- as.matrix(stemp5)
		stemp7 <- as.matrix(stemp7)
        sigmean <- sigmean + stemp5								#	diagonal of the VAR covariance matrix
        cormean <- cormean + stemp7								#	off-diagonal elements of the VAR cov matrix
        sig2mo <- sig2mo + stemp5^2
        cor2mo <- cor2mo + stemp7^2

	  } # -- End of "saving after burn-in results" 		
	
  } # - End of "for 1" main Gibbs loop (for irep = 1:N_mcmc+N_burn)
print('Time-processing elapsed:')
proc.time() - tic
#	--------------------------	End of MCMC sampling	-----------------------------
#	---------------------------------------------------------------------------------
  
#	Output
output <- list(Bt_postmean,At_postmean,Sigt_postmean,Qmean,Smean,Wmean,sigmean,cormean)
return(output)
}
#
# -- End of "tvpvar.R" -- 
