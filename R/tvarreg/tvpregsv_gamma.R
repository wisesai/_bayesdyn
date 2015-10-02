tvpregsv_gamma <- function(y,x,w,beta_sim,alpha_sim,h_sim,gamma_sim,a_gamma,b_gamma){
#
#	Inputs:
#   Observables
#	    y    		- response variable (time-series vector of dimension T x 1)
#		x			- covariates (time-series matrix of dimension T x k)
#		w			- covariates (time-series matrix of dimension T x p)
#		beta_sim	- current state for the 'beta' parameter
#		alpha_sim	- current state for the '(alpha_t)' hidden variables
#		h_sim		- current state for the '(h-t)' hidden variables
#		a_gamma		- hyper shape parameter of the prior distribution for 'gamma' 
#		b_gamma		- hyper scale parameter of the prior distribution for 'gamma' 
#
#   Author: 	Juan Carlos Martinez-Ovando 
#				(Email: juan.martinez@banxico.org.mx)
#				(Email: JC.Martinez.Ovando@gmail.com)
#
#	Revision:		May 17, 2015
#

#	Repository
gamma_sim_new <- matrix(NaN, 1, 1)
T <- dim(y)[1]

#	Defining temporal variables
yxw <- y - x%*%beta_sim - rowprodmatrix(w,t(alpha_sim))
Sigmat_inv <- diag(c(1/exp(h_sim)))

#	Updating parameters
a_gamma_post <- a_gamma + T/2
b_gamma_post <- b_gamma + 0.5 * t(yxw)%*%Sigmat_inv%*%yxw
gamma_sim_new <- as.matrix(rgamma(n=1, shape=a_gamma_post, scale=1/b_gamma_post))
#	Note: E(gamma) = a_gamma_post / b_gamma_post 
#gamma_sim_new <- 1/gamma_sim_new

#	Output
return(gamma_sim_new)
}
#
#   --  End of "tvpregsv_gamma.R" --