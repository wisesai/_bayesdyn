tvpvar_prior <- function(prior_type,Y,N_training_prior,M,p_lags){
#
#	tvpvar_primiceri.R
#	Time Series prior under Primiceri's approach
#
#	Original code was written by D. Korobilis and adpated by J.C. Martinez-Ovando (2015)
#
#	Author: 	Juan~Carlos Martínez-Ovando
#	Email:		juan.martinez@banxico.org.mx
#				jc.martinez.ovando@gmail.com
#

if(prior_type == 1){
	#	To set up training sample prior a-la Primiceri, use the following subroutine
	ts_prior_out <- tvpvar_primiceri(Y,N_training_prior,M,p_lags)
	B_OLS <- ts_prior_out[[1]]
	VB_OLS <- ts_prior_out[[2]]
	A_OLS <- ts_prior_out[[3]]
	VA_OLS <- ts_prior_out[[5]]
	sigma_OLS <- ts_prior_out[[4]]
  }else if(prior_type == 1){
	#	----	Or use uninformative values	----
	B_OLS <- matrix(0, K, 1)
	VB_OLS <- diag(K)
	A_OLS <- matrix(0, numa, 1)
	VA_OLS <- diag(numa)
	sigma_OLS <- 0*matrix(1, M, 1)
  }else{
	stop('tvpvar: Error in the specifiction of the type for the prior distribution.')
  }
  
#	Output
output <- list(B_OLS,VB_OLS,A_OLS,VA_OLS,sigma_OLS)
return(output) 
}
#
#  --  End of "tvpvar_prior.R"  --