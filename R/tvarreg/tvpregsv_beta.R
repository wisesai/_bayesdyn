tvpregsv_beta <- function(beta_sim,y,x,w,alpha_sim,gamma_sim,h_sim,m_beta,C_beta){
#
#	Inputs:
#   Observables
#		beta_sim	- current state for the 'beta' parameter
#	    y    		- response variable (time-series vector of dimension T x 1)
#		x			- covariates (time-series matrix of dimension T x k)
#		w			- covariates (time-series matrix of dimension T x p)
#		gamma_sim	- current state for the 'gamma' parameter
#		alpha_sim	- current state for the hidden '(alpha_t)' parameters
#		h_sim		- current state for the hidden parameters '(h_t)'
#   Hyper-priors
#		m_beta 		- mean parameter of the prior distribution for 'beta'
#		C_beta		- covariance matrix of the prior distribution for 'beta'
#
#   Author: 	Juan Carlos Martinez-Ovando 
#				(Email: juan.martinez@banxico.org.mx)
#				(Email: JC.Martinez.Ovando@gmail.com)
#
#	Revision:		May 17, 2015
#

#	Repository
k <- dim(beta_sim)[1]
beta_sim_new <- matrix(NaN, k, 1)

#	Scanning
yw <- y - rowprodmatrix(w,t(alpha_sim))
Sigmmat_inv <- diag(c(1/(as.numeric(gamma_sim) * exp(h_sim))))
C_beta_inv <- solve(C_beta)

#	Posterior
C_post_inv <- t(x)%*%Sigmmat_inv%*%x + C_beta_inv
C_post <- solve(C_post_inv)
m_post <- C_post %*% ( t(x)%*%Sigmmat_inv%*%yw - C_beta_inv%*%m_beta )
beta_sim_new <- as.matrix(mvrnorm( n=1, mu=m_post, Sigma=C_post))

#	Output
return(beta_sim_new)
}
#
#   --  End of "tvpregsv_beta.R" --