tvpregsv_eta <- function(eta_sim, h_sim,z,h_0,sigmah_sim,varphi_sim,phi_sim,mu_eta,var_eta,N_slice){
#
#	Inputs:
#   Observables
#		eta_sim		- current state of the parameter 'eta'
#		h_sim		- current state for the hidden '(h_t)' parameters
#		z			- covariates (time-series matrix of dimension T x q)
#		h_0			- initial value of the hidden markov '(h_t)'
#		sigmah_sim	- current state for scale parameter for the hidden markov '(h_t)'
#		varphi_sim	- current state for the parameter 'varphi'
#		phi_sim		- current state for the parameter 'phi'
#		mu_eta		- mean parameter for the prior for 'eta'
#		var_eta		- variance parameter for the prior for 'eta'
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
eta_sim_new <- NaN * eta_sim
T <- dim(h_sim)[1]

#	Defining temporal variables
h_sim_varphi <- NaN * h_sim
h_sim_varphi[1,] <- varphi_sim * h_0
h_sim_varphi[2:T,] <- varphi_sim * h_sim[1:(T-1),] 
hz_temp <- (h_sim - z%*%phi_sim - h_sim_varphi) / as.numeric(1-varphi_sim) 
sigmah_slice <- as.numeric(sigmah_sim / as.numeric(1-varphi_sim)^2)

#	Slice sampler steps
eta_ini <- as.numeric(eta_sim)
uni.slice.calls <<- 0	# Number of calls of the slice sampling function
uni.slice.evals <<- 0	# Number of density evaluations done in these calls
eta_sim_new <- uni.slice(eta_ini, function (eta){ -1/(2*sigmah_slice)*sum((hz_temp - eta)^2) 
												  -1/(2*var_eta)*(eta-mu_eta)^2 }, m=13, lower=-1e10, upper=1e10)

#	Output
return(eta_sim_new)
}
#
#   --  End of "tvpregsv_eta.R" --