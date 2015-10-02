#
#	Demo for "bvar.R"
#

rm(list=ls())

#	Paths
path.code <- 'L:/JCMO.Research/Coding/_r/_my_packages/_bayesdyn/bayesdyn/R'
path.data <- 'L:/JCMO.Research/Coding/_r/_my_packages/_bayesdyn/bayesdyn/R/bvar'
path.out <- 'L:/JCMO.Research/Coding/_r/_my_packages/_bayesdyn/bayesdyn/R/bvar'

#	Packages
#install.packages('mvtnorm')
library('mvtnorm')
#install.packages('pracma')
library('pracma')

#	Source code
source(paste(path.code,"/utilities/finding.R",sep=""))
source(paste(path.code,"/utilities/indexownlags.R",sep=""))
source(paste(path.code,"/utilities/matrixcol.R",sep=""))
source(paste(path.code,"/utilities/mlag2.R",sep=""))
source(paste(path.code,"/utilities/reshaper.R",sep=""))
source(paste(path.code,"/utilities/reshapes.R",sep=""))
source(paste(path.code,"/utilities/squeeze.mean.R",sep=""))
source(paste(path.code,"/utilities/wishartrnd.R",sep=""))

source(paste(path.code,"/bvar/bvar.R",sep=""))

set.seed(2)

#	------------------------------  Data  -----------------------------------
#	Load Quarterly US data on inflation, unemployment and interest rate, 
#	1953:Q1 - 2006:Q3
Yraw <- read.csv(paste(path.data,"/ydata.csv",sep=""), header = TRUE, row.names = 1)
Yraw <- as.matrix(Yraw)
#Yrawts <- ts(Yraw, frequency = 4, start = c(1953,1))
#plot(Yrawts)

#	Note:	'Yraw' is a matrix with T rows by M columns, where T is the number 
#			of time series observations (usually months or quarters), 
#			while M is the number of VAR dependent macro variables.

#	--------------------------  Preliminaries  -------------------------------
#	Define specification of the VAR model
constant <- 1				#	1: if you desire intercepts, 0: otherwise 
p_lags <- 2					#	Number of lags on dep_lagsendent variables

forecasting <- 1			#	1: Compute h-step ahead predictions, 0: no prediction

forecast_method <- 1		#	0: Direct forecasts, 1: Iterated forecasts

N_pred_thin <- 10			#	Number of times to obtain a draw from the predictive 
							#	density, for each generated draw of the parameters                     

N_pred_h <- 13				#	Number of forecast periods

impulses <- 1				#	1: compute impulse responses, 0: no impulse responses

impulses_h <- 21			#	Horizon to compute impulse responses

#	Set prior for BVAR model:
prior <- 1					#	prior = 1 --> Indepependent Normal-Whishart Prior
							#	prior = 2 --> Indepependent Minnesota-Whishart Prior

#	Gibbs-related preliminaries
N_mcmc <- 100				#	Final number of draws to save
N_burn <- 20				#	Draws to discard (burn-in)
N_print <- 10				#	Print on the screen every "N_print"-th iteration

#-------------------------------------------------------------------------------------
#									Gibbs sampler
#-------------------------------------------------------------------------------------
bvar_out <- bvar(Yraw, constant, p_lags, forecasting, forecast_method, N_pred_h, impulses, impulses_h, prior, N_mcmc, N_burn, N_print)
alpha_draws <- bvar_out[[1]]
ALPHA_draws <- bvar_out[[2]]
SIGMA_draws <- bvar_out[[3]]

#	Posterior mean of parameters:
ALPHA_mean <- squeeze.mean(ALPHA_draws)									#	posterior mean of ALPHA
SIGMA_mean <- squeeze.mean(SIGMA_draws)									#	posterior mean of SIGMA

#	mean prediction and log predictive likelihood
if(forecasting == 1){
    Y_pred_mean <- colMeans(Y_pred)
    log_PL <- mean((log(PL)))
	#	This are the true value of Y at T+N_pred_h:
	if(forecast_method == 0){
		true_value <- Y1[T+1, ]
	  }else if( forecast_method == 1 ) {
		true_value <- Y1[T+N_pred_h, ]
	  }
	#	(subsequently you can easily also get MSFE and MAFE)
	}

#	You can also get other quantities, like impulse responses
if(impulses == 1){
    qus <- c(0.1, .5, .90)
    imp_resp <- squeeze(quantile(imp,qus));
  }

#
# -- End of "BVAR_GIBBS.R" --