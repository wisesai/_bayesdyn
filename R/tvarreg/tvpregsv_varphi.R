tvpregsv_varphi <- function(varphi_sim, h_sim,z,h_0,sigmah_sim,eta_sim,phi_sim,a_varphi,b_varphi,N_slice){
#
#	Inputs:
#   Observables
#		varphi_sim	- current state for the parameter 'varphi'
#		h_sim		- current state for the hidden '(h_t)' parameters
#		z			- covariates (time-series matrix of dimension T x q)
#		h_0			- initial value of the hidden markov '(h_t)'
#		sigmah_sim	- current state for scale parameter for the hidden markov '(h_t)'
#		eta_sim		- current state of the parameter 'eta'
#		phi_sim		- current state for the parameter 'phi'
#		a_varphi	- mean parameter for the prior for 'varphi'
#		b_varphi	- variance parameter for the prior for 'varphi'
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

#	Repository for new simulations
varphi_sim_new <- NaN * varphi_sim
T <- dim(h_sim)[1] 

#	Defining temporal variables
hetaz <- h_sim - as.numeric(eta_sim) - z%*%phi_sim
heta <- rbind(as.matrix(h_0),as.matrix(h_sim[1:(T-1),])) - as.numeric(eta_sim)
Sigmah_inv <- (1/as.numeric(sigmah_sim)) * diag(T)

#	Slice sampler steps
varphi_ini <- as.numeric(varphi_sim)
varphi_sim_new <- uni.slice(varphi_ini, function (varphi){ -0.5*(t(hetaz - varphi*heta)%*%Sigmah_inv%*%(hetaz - varphi_ini*heta))
														   +(a_varphi-1)*log((varphi+1)/2) 
														   +(b_varphi-1)*log(1-(varphi+1)/2) }, m=13, lower=(-1+1e-10), upper=(1-1e-10))

#	Output
return(varphi_sim_new)
}
#
#   --  End of "tvpregsv_varphi.R" --