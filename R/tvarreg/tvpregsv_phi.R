tvpregsv_phi <- function(phi_sim,h_sim,z,h_0,sigmah_sim,eta_sim,varphi_sim,m_phi,C_phi){
#
#	Inputs:
#   Observables
#		phi_sim		- current state for the hidden 'phi' parameter
#		h_sim		- current state for the hidden '(h_t)' parameters
#		z			- covariates (time-series matrix of dimension T x q)
#		h_0			- initial value of the hidden markov '(h_t)'
#		sigmah_sim	- current state for scale parameter for the hidden markov '(h_t)'
#		eta_sim		- current state of the parameter 'eta'
#		varphi_sim	- current state for the parameter 'varphi'
#		m_phi		- mean parameter for the prior for 'phi'
#		C_phi		- covariance parameter for the prior for 'phi'
#
#   Author: 	Juan Carlos Martinez-Ovando 
#				(Email: juan.martinez@banxico.org.mx)
#				(Email: JC.Martinez.Ovando@gmail.com)
#
#	Revision:		May 17, 2015
#

#	Repository for new simulations
phi_sim_new <- NaN * phi_sim
T <- dim(h_sim)[1] 

#	Defining temporal variables
heta <- h_sim - as.numeric(eta_sim) - as.numeric(varphi_sim)*(rbind(as.matrix(h_0),as.matrix(h_sim[1:(T-1),])) - as.numeric(eta_sim))
Sigmah_inv <- (1 / as.numeric(sigmah_sim)) * diag(T)
C_phi_inv <- solve(C_phi)

#	Posterior
C_post_inv <- t(z)%*%Sigmah_inv%*%z + C_phi_inv
C_post <- solve(C_post_inv)
m_post <- C_post %*% ( t(z)%*%Sigmah_inv%*%heta - C_phi_inv%*%m_phi )
phi_sim_new <- as.matrix(mvrnorm( n=1, mu=m_post, Sigma=C_post))

#	Output
return(phi_sim_new)
}
#
#   --  End of "tvpregsv_phi.R" --