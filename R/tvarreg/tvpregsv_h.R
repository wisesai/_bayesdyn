tvpregsv_h <- function(h_sim,y,x,w,z,h_0,sigmah_sim,beta_sim,alpha_sim,gamma_sim,eta_sim,varphi_sim,phi_sim,N_slice){
#
#	Inputs:
#   Observables
#		h_sim		- current state for the hidden '(h_t)' parameters
#	    y    		- response variable (time-series vector of dimension T x 1)
#		x			- covariates (time-series matrix of dimension T x k)
#		w			- covariates (time-series matrix of dimension T x p)
#		z			- covariates (time-series matrix of dimension T x q)
#		h_0			- initial value of the hidden markov '(h_t)'
#		sigmah_sim	- current state for scale parameter for the hidden markov '(h_t)'
#		beta_sim	- current state for the 'beta' parameter
#		alpha_sim	- current state for the '(alpha_t)' hidden variables
#		gamma_sim	- current state for the 'gamma' parameter
#		eta_sim		- current state for the parameter 'eta'
#		varphi_sim	- current state for the parameter 'varphi'
#		phi_sim		- current state for the parameter 'phi'
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
h_sim_new <- NaN * h_sim
T <- dim(y)[1]

#	Defining temporal variables
yxw <- y - x%*%beta_sim - rowprodmatrix(w,t(alpha_sim)) 
hetaz <- h_sim - as.numeric(eta_sim) - as.numeric(varphi_sim)*(rbind(as.matrix(h_0),as.matrix(h_sim[1:(T-1),])) - as.numeric(eta_sim)) - z%*%phi_sim

#	Markovian recursion
t <- 1
for(t in 1:1){
	sigmah_slice <- as.numeric(1/(2*sigmah_sim))
	#	Part 1: Observation equation
	yxw2_part1 <- (yxw[t,]^2)/(2*gamma_sim)
	#	Part 2: Initial for h_{t-1}
	mu_phi_part2 <- eta_sim + as.numeric(varphi_sim*(h_0-eta_sim)) + as.numeric(z[t, ]%*%phi_sim)
	#	Part 3: Initial for h_{t+1}
	mu_phi_part3 <- (h_sim[(t+1),] - eta_sim) - as.numeric(varphi_sim*eta_sim) + as.numeric(z[t, ]%*%phi_sim)
	#	Slice sampler step
	h_ini <- as.numeric(h_sim[t, ])
	uni.slice.calls <<- 0	# Number of calls of the slice sampling function
	uni.slice.evals <<- 0	# Number of density evaluations done in these calls
	h_sim_new[t, ] <- uni.slice(h_ini, function (h){  -0.5 * h
													  -yxw2_part1 / exp(h)
													  -sigmah_slice * (h - mu_phi_part2)^2 
													  -sigmah_slice * (mu_phi_part3 - varphi_sim*h)^2 }, m=N_slice, lower=-1e10, upper=1e10)
  }

t <- 2
for(t in 2:(T-1)){
	sigmah_slice <- as.numeric(1/(2*sigmah_sim))
	#	Part 1: Observation equation
	yxw2_part1 <- (yxw[t,]^2)/(2*gamma_sim)
	#	Part 2: Initial for h_{t-1}
	mu_phi_part2 <- eta_sim + as.numeric(varphi_sim*(h_sim_new[(t-1),]-eta_sim)) + as.numeric(z[t, ]%*%phi_sim)
	#	Part 3: Initial for h_{t+1}
	mu_phi_part3 <- (h_sim[(t+1),] - eta_sim) - as.numeric(varphi_sim*eta_sim) + as.numeric(z[t, ]%*%phi_sim)
	#	Slice sampler step
	h_ini <- as.numeric(h_sim[t, ])
	uni.slice.calls <<- 0	# Number of calls of the slice sampling function
	uni.slice.evals <<- 0	# Number of density evaluations done in these calls
	h_sim_new[t, ] <- uni.slice(h_ini, function (h){  -0.5 * h
													  -yxw2_part1 / exp(h)
													  -sigmah_slice * (h - mu_phi_part2)^2 
													  -sigmah_slice * (mu_phi_part3 - varphi_sim*h)^2 }, m=N_slice, lower=-1e10, upper=1e10)
  }

t <- T
for(t in T:T){
	sigmah_slice <- as.numeric(1/(2*sigmah_sim))
	#	Part 1: Observation equation
	yxw2_part1 <- (yxw[t,]^2)/(2*gamma_sim)
	#	Part 2: Initial for h_{t-1}
	mu_phi_part2 <- eta_sim + as.numeric(varphi_sim*(h_sim_new[(t-1),]-eta_sim)) + as.numeric(z[t, ]%*%phi_sim)
	#	Slice sampler step
	h_ini <- as.numeric(h_sim[t, ])
	uni.slice.calls <<- 0	# Number of calls of the slice sampling function
	uni.slice.evals <<- 0	# Number of density evaluations done in these calls
	h_sim_new[t, ] <- uni.slice(h_ini, function (h){  -0.5 * h
													  -yxw2_part1 / exp(h)
													  -sigmah_slice * (h - mu_phi_part2)^2 }, m=N_slice, lower=-1e10, upper=1e10)
  }

#	Output
return(h_sim_new)
}
#
#   --  End of "tvpregsv_h.R" --