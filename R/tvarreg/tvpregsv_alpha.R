tvpregsv_alpha <- function(alpha_sim,y,x,w,alpha_0,beta_sim,gamma_sim, h_sim,sigmaa_sim){
#
#	Inputs:
#   Observables
#		alpha_sim	- current state for the hidden '(alpha_t)' parameters
#	    y    		- response variable (time-series vector of dimension T x 1)
#		x			- covariates (time-series matrix of dimension T x k)
#		w			- covariates (time-series matrix of dimension T x p)
#		alpha_0		- initial value of the hidden markov '(alpha_t)'
#		beta_sim	- current state for the 'beta' parameter
#		gamma_sim	- current state for the 'gamma' parameter
#		h_sim		- current state for the hidden parameters '(h_t)'
#		sigmaa_sim	- current state for the 'sigma_a' parameter (variance of the evolution for 'alpha_t')
#
#   System
#	    N_mcmc    	- length of the Gibbs sampler
#	    N_slice   	- length of the slice sampler (typically 1)
#
#   Author: 	Juan Carlos Martinez-Ovando 
#				(Email: juan.martinez@banxico.org.mx)
#				(Email: JC.Martinez.Ovando@gmail.com)
#
#	Revision:		May 17, 2015
#

#	Repository for new simulations
alpha_sim_new <- NaN * as.matrix(alpha_sim)
T <- dim(y)[1]

#	Defining temporal variables
yx <- y - x%*%beta_sim
Sigmat_inv <- diag(c(1/(as.numeric(gamma_sim) * exp(h_sim))))
Sigmaa_inv <- solve(sigmaa_sim)

#	Markovian recursion
t <- 1
for(t in 1:1){
	C_updated_inv <- as.matrix(w[t, ])%*%Sigmat_inv[t,t]%*%t(w[t, ]) + Sigmaa_inv
	C_updated <- solve(C_updated_inv)
	mu_updated <- C_updated %*% ( as.matrix(w[t, ])%*%Sigmat_inv[t,t]%*%t(yx[t,]) + Sigmaa_inv%*%alpha_0 )
	alpha_sim_new[ ,t] <- as.matrix(mvrnorm(n = 1, mu=mu_updated, Sigma=C_updated))
  }

t <- 2
for(t in 2:(T-1)){
	C_updated_inv <- as.matrix(w[t, ])%*%Sigmat_inv[t,t]%*%t(w[t, ]) + Sigmaa_inv + Sigmaa_inv
	C_updated <- solve(C_updated_inv)
	mu_updated <- C_updated %*% ( as.matrix(w[t, ])%*%Sigmat_inv[t,t]%*%t(yx[t,]) + Sigmaa_inv%*%alpha_sim_new[,(t-1)] + Sigmaa_inv%*%alpha_sim[,(t+1)] )
	alpha_sim_new[ ,t] <- as.matrix(mvrnorm(n = 1, mu=mu_updated, Sigma=C_updated))
  }

t <- T
for(t in T:T){
	C_updated_inv <- as.matrix(w[t, ])%*%Sigmat_inv[t,t]%*%t(w[t, ]) + Sigmaa_inv
	C_updated <- solve(C_updated_inv)
	mu_updated <- C_updated %*% ( as.matrix(w[t, ])%*%Sigmat_inv[t,t]%*%t(yx[t,]) + solve(sigmaa_sim)%*%alpha_sim_new[,(t-1)] )
	alpha_sim_new[ ,t] <- as.matrix(mvrnorm(n = 1, mu=mu_updated, Sigma=C_updated))
  }

#	Output
return(alpha_sim_new)
}
#
#   --  Eend of "tvpregsv_alpha.R" --