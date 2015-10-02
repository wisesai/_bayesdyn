tvpregsv <- function(y,x,w,z,N_mcmc,N_slice){
#
#   This function computes a Gibbs sampler of a time-varying parameter regression model with
#	dependent stochastic volatility component
#
#	Inputs:
#   Observables
#	    y    		- response variable (time-series vector of dimension T x 1)
#		x			- covariates (time-series matrix of dimension T x k)
#		w			- covariates (time-series matrix of dimension T x p)
#		z			- covariates (time-series matrix of dimension T x q)
#   System
#	    N_mcmc    	- length of the Gibbs sampler
#	    N_slice   	- length of the slice sampler (typically 1)
#
#   Outputs:
#		beta_mcmc	- regression parameters (vector of dimension 2 x N_mcmc)			
#		alpha_mcmc	- time-varying effects on mean (vector of dimension T x N_mcmc)			
#		gamma_mcmc	- volatility parameter (vector of dimension 1 x N_mcmc)			
#		h_mcmc		- time-varying volatilities (vector of dimension T x N_mcmc)			
#		eta_mcmc	- time-varying volatilities (vector of dimension 1 x N_mcmc for the 'eta' parameter)			
#		varphi_mcmc	- stochastic volatility parameters (vector of dimension 1 x N_mcmc for the 'varphi' parameter)			
#		phi_mcmc	- stochastic volatility parameters (vector of dimension q x N_mcmc for the 'phi' parameter)			
#		sigmah_mcmc	- variance parameter for mean (vector of dimension 1 x N_mcmc)			
#
#   Author: 	Juan Carlos Martinez-Ovando 
#				(Email: juan.martinez@banxico.org.mx)
#				(Email: JC.Martinez.Ovando@gmail.com)
#
#	Revision:		May 17, 2015
#

#	Initial scanning 
T <- dim(y)[1]
k <- dim(x)[2]
p <- dim(w)[2]
q <- dim(z)[2]

#	Setting up storage arrays
beta_mcmc <- NaN * matrix(1, k, N_mcmc)				#	regression parameter
alpha_mcmc <- NaN * array(1, c(p,T,N_mcmc))			#	hidden time-varying regresion variable
gamma_mcmc <- NaN * matrix(1, 1, N_mcmc)			#	hidden time-varying volatility
h_mcmc <- NaN * matrix(1, T, N_mcmc)				#	hidden time-varying volatility
sigmah_mcmc <- NaN * matrix(1, 1, N_mcmc)			#	volatility parameter for hidden volatility variable
eta_mcmc <- NaN * matrix(1, 1, N_mcmc)				#	volatility - 'eta' parameter
varphi_mcmc <- NaN * matrix(1, 1, N_mcmc)			#	volatility - AR parameter
phi_mcmc <- NaN * matrix(1, q, N_mcmc)				#	volatility - regression parameter
sigmaa_mcmc <- NaN * matrix(1, 1,N_mcmc)			#	volatility parameter for hidden time-varying regression variable

#	-------------------------------------------------------
#	-------------------------------------------------------
#		Setting up initial values for hybrid-Gibbs sampler
#		NOTE: This are to be empirically defined...
#	-------------------------------------------------------
#	-------------------------------------------------------
if(type_ini == 1){
	#	Informative initial values
	x_temp <- cbind(x, w)
	betaalpha_ini <- solve(t(x_temp) %*% x_temp) %*% t(x_temp) %*% y
	betaalpha_ini
	beta_sim <- as.matrix(betaalpha_ini[1:k, ])							#	regression parameter
	alpha_sim <- repmat(as.matrix(betaalpha_ini[(k+1):(k+p), ]),1,T) 	#	hidden time-varying regresion variable
	sigmaa_sim <- 0.3 * diag(p)												#	volatility parameter for hidden time-varying regression variable
	error_temp <- y - x_temp%*%betaalpha_ini
	gamma_sim <- 1.0
	error_var <- NaN * error_temp
	t <- 1
	for(t in 1:T){
		if(t < 9){
			error_var[t] <- var(error_temp[t:(t+9)])
		  }else if( (t >= 9) & (t <= (T-9)) ){
			error_var[t] <- var(error_temp[(t-9):(t+9)])
		  }else{
			error_var[t] <- var(error_temp[(t-9):t])
		  }
	  }
	h_sim <- log(error_var/gamma_sim)									#	hidden time-varying volatility
	h_sim_lag <- rbind(matrix(1,1,1), as.matrix(h_sim[2:T,]))
	plot(h_sim[2:T],h_sim[1:(T-1)])
	sigmah_sim <- matrix(var(h_sim-mean(h_sim)), 1, 1)					#	volatility parameter for hidden volatility variable
	eta_sim <- mean(h_sim)
	h_temp <- cbind( (h_sim_lag-eta_sim-eta_sim), z)
	colnames(h_temp) <- c("h_t_lag","z_t")
	varphiphi_ini <- solve(t(h_temp) %*% h_temp) %*% t(h_temp) %*% error_temp
	eta_sim <- as.matrix(eta_sim)										#	volatility - 'eta' parameter
	varphi_sim <- as.matrix(varphiphi_ini[1, ])							#	volatility - 'varphi' parameter
	phi_sim <- as.matrix(varphiphi_ini[(1+1):(1+q),])					#	volatility - 'phi' parameter
  }else if(type_ini == 0){
	#	Vague initial values
	beta_sim <- matrix(0,k,1)											#	regression parameter
	alpha_sim <- repmat(matrix(0,p,1),1,T) 								#	hidden time-varying regresion variable
	sigmaa_sim <- 0.3*diag(p)												#	volatility parameter for hidden time-varying regression variable
	gamma_sim <- 1.0
	h_sim <- matrix(0,T,1)												#	hidden time-varying volatility
	sigmah_sim <- matrix(1, 1, 1)										#	volatility parameter for hidden volatility variable
	eta_sim <- mean(h_sim)
	eta_sim <- as.matrix(eta_sim)										#	volatility - 'eta' parameter
	varphi_sim <- matrix(0, 1, 1)										#	volatility - 'varphi' parameter
	phi_sim <- matrix(0, q, 1)											#	volatility - 'phi' parameter
  }
  #
  #		Hidden variables
  #
  alpha_0 <- as.matrix(rowMeans(alpha_sim))								#	initial value of the Markov chain for the hidden variables '(alpha_t)'
  h_0 <- mean(h_sim)													#	initial value of the Markov chain for the hidden variables '(h_t)'

#	-------------------------------------------------------
#	-------------------------------------------------------
#		Prior specification
#	-------------------------------------------------------
#	-------------------------------------------------------
m_beta <- matrix(0,k,1)									#	mean parameter for the 'beta' parameter
C_beta <- diag(k)										#	covariance matrix for the 'beta' parameter

a_gamma <- 0.3											#	hyper shape-parameter of prior distribution for 'gamma' parameter
b_gamma <- 0.3											#	hyper scale-parameter of prior distribution for 'gamma' parameter

mu_eta <- 0												#	mean parameter for the prior for 'eta' (which is Gaussian)
var_eta <- 1											#	variance parameter for the prior for 'eta' (which is Gaussian)

a_varphi <- 3											#	shape 1 parameter for the prior for 'varphi' (which is Generalized beta at (-1,1) )
b_varphi <- 3											#	shape 2 parameter for the prior for 'varphi' (which is Generalized beta at (-1,1) )

m_phi <- matrix(0, q,1)									#	mean parameter for the prior for 'phi' (which is Gaussian )
C_phi <- diag(q)										#	covariance parameter for the prior for 'phi' (which is Gaussian )

a_sigmah <- 0.3											#	shape-parameter of prior distribution for 'sigmah' parameter
b_sigmah <- 0.3											#	scale-parameter of prior distribution for 'sigmah' parameter

#---------------------------------------------------	
#---------------------------------------------------	
#				MCMC Sampling
#---------------------------------------------------	
#---------------------------------------------------	
print('=== Hybrid Gibbs sampler for tvpregsv.r ===');
it_mcmc = 1
for(it_mcmc in 1:N_mcmc){
	print( paste("MCMC: ", it_mcmc, " of ", N_mcmc, sep=""))

	#	--------------------------
	#		Sampling parameters
	#	--------------------------
    # #
	# # Sampling 'beta'
	# #
	beta_mcmc[ ,it_mcmc] <- beta_sim					#	parameters
	beta_sim <- tvpregsv_beta(beta_sim,y,x,w,alpha_sim,gamma_sim,h_sim,m_beta,C_beta)

	# #
	# # Sampling 'gamma'
	# #
	gamma_mcmc[ , it_mcmc] <- gamma_sim
	gamma_sim <- tvpregsv_gamma(y,x,w,beta_sim,alpha_sim,h_sim,gamma_sim,a_gamma,b_gamma)

	# #
	# # Sampling 'eta'
	# #
	eta_mcmc[ ,it_mcmc] <- eta_sim
	eta_sim <- tvpregsv_eta(eta_sim, h_sim,z,h_0,sigmah_sim,varphi_sim,phi_sim,mu_eta,var_eta,N_slice)

	# #
	# # Sampling 'varphi'
	# #
	varphi_mcmc[ ,it_mcmc] <- varphi_sim
	# varphi_sim <- tvpregsv_varphi(varphi_sim, h_sim,z,h_0,sigmah_sim,eta_sim,phi_sim,a_varphi,b_varphi,N_slice)

	# #
	# # Sampling 'phi'
	# #
	phi_mcmc[ ,it_mcmc] <- phi_sim
	phi_sim <- tvpregsv_phi(phi_sim,h_sim,z,h_0,sigmah_sim,eta_sim,varphi_sim,m_phi,C_phi)

	# #
	# # Sampling 'sigmah'
	# #
	sigmah_mcmc[ ,it_mcmc] <- sigmah_sim
	sigmah_sim <- tvpregsv_sigmah(h_sim,z,h_0,eta_sim,varphi_sim,phi_sim,sigmah_sim,a_sigmah,b_sigmah)

	#	---------------------------------
	#		Sampling hidden variables
	#	---------------------------------
	# #
	# # Sampling 'alpha'
	# #
	alpha_mcmc[ , , it_mcmc] <- alpha_sim
	alpha_sim <- tvpregsv_alpha(alpha_sim,y,x,w,alpha_0,beta_sim,gamma_sim, h_sim,sigmaa_sim)
	
	# #
	# # Sampling 'h'
	# #
	h_mcmc[ ,it_mcmc] <- h_sim
	# h_sim <- tvpregsv_h(h_sim,y,x,w,z,h_0,sigmah_sim,beta_sim,alpha_sim,gamma_sim,eta_sim,varphi_sim,phi_sim,N_slice)
		
  }

#	Output
output <- list(beta_mcmc,alpha_mcmc,gamma_mcmc,h_mcmc,eta_mcmc,varphi_mcmc,phi_mcmc,sigmah_mcmc)
return(output)
}
#
#   --  End of "tvpregsv.R" --