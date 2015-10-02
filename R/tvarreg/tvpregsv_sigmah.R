tvpregsv_sigmah <- function(h_sim,z,h_0,eta_sim,varphi_sim,phi_sim,sigmah_sim,a_sigmah,b_sigmah){
#
#	Inputs:
#   Observables
#		h_sim		- current state for the hidden '(h_t)' parameters
#		z			- covariates (time-series matrix of dimension T x q)
#		h_0			- initial value of the hidden markov '(h_t)'
#		eta_sim		- current state for the parameter 'eta'
#		varphi_sim	- current state for the parameter 'varphi'
#		phi_sim		- current state for the parameter 'phi'
#		sigmah_sim	- current state for the parameter 'sigmah'
#		a_sigmah	- (hyperparameter) shape parameter of the prior for 'sigmah'
#		b_sigmah	- (hyperparameter) scale parameter of the prior for 'sigmah'
#
#   System
#	    N_slice   	- length of the slice sampler (typically 1)
#
#   Author: 	Juan Carlos Martinez-Ovando 
#				(Email: juan.martinez@banxico.org.mx)
#				(Email: JC.Martinez.Ovando@gmail.com)
#
#	Revision:		May 17, 2015
#

#	Repository
sigmah_sim_new <- NaN * sigmah_sim
T <- dim(h_sim)[1] 

#	Defining temporal variables
hetaz <- h_sim - as.numeric(eta_sim) - as.numeric(varphi_sim)*(rbind(as.matrix(h_0),as.matrix(h_sim[1:(T-1),])) - as.numeric(eta_sim)) - z%*%phi_sim

#	Updating parameters
a_sigmah_post <- a_sigmah + T/2
b_sigmah_post <- b_sigmah + 0.5 * t(hetaz)%*%hetaz
sigmah_sim_new <- as.matrix(rgamma(n=1, shape=a_sigmah_post, scale=1/b_sigmah_post))
#	Note: E(gamma) = a_sigmah_post / b_sigmah_post 
#sigmah_sim_new <- 1/sigmah_sim_new

#	Output
return(sigmah_sim_new)
}
#
#   --  End of "tvpregsv_sigmah.R" --