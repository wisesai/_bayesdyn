#
#	Demo for "dynfactor.R"
#

rm(list=ls())

#	Paths
path.code <- 'L:/JCMO.Research/Coding/_r/_my_packages/_bayesdyn/bayesdyn/R/dynfactor'
path.data <- 'L:/JCMO.Research/Coding/_r/_my_packages/_bayesdyn/bayesdyn/R/dynfactor'
path.out <- 'L:/JCMO.Research/Coding/_r/_my_packages/_bayesdyn/bayesdyn/R/dynfactor'

#	Packages
#install.packages('mvtnorm')
library('mvtnorm')
#install.packages('pracma')
library('pracma')
#install.packages('MASS')
library('MASS')

#	Source code
source(paste(path.code,"/dfmsim.R",sep=""))
source(paste(path.code,"/mlag.R",sep=""))
source(paste(path.code,"/pcext.R",sep=""))
source(paste(path.code,"/reshaper.R",sep=""))
source(paste(path.code,"/reshapes.R",sep=""))
source(paste(path.code,"/varfreq.R",sep=""))
source(paste(path.code,"/squeeze.mean.R",sep=""))
source(paste(path.code,"/dynfactor.R",sep=""))

set.seed(2)

#	---------------------------LOAD DATA--------------------------------------
#		Simulated dataset
X <- dfmsim()
demean <- 1

#	Number of factors & lags in B(L):
K <- 3
lags <- 1

#	----------------------------PRELIMINARIES---------------------------------
#	Set some Gibbs - related preliminaries
N_mcmc <- 100							#	Number of draws
N_burn <- 0							#	Number of burn-in-draws
N_thin <- 1							#	Consider every N_thin-th draw (N_thin value)
N_print <- 1						#	Number of jumps for printing out iterations

#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#									Gibbs sampler
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
dynfactor_out <- dynfactor(X, demean, K, lags, N_mcmc, N_burn, N_thin, N_print)
Ldraw <- dynfactor_out[[1]]
Bdraw <- dynfactor_out[[2]]
Qdraw <- dynfactor_out[[3]]
Fdraw <- dynfactor_out[[4]]
Rdraw <- dynfactor_out[[5]]

#	Do thining in case of high correlation
thin_val <- seq(1,((N_mcmc-N_burn)/N_thin), by = N_thin)
Ldraw <- Ldraw[thin_val, , ]
Bdraw <- Bdraw[thin_val, , , ]
Qdraw <- Qdraw[thin_val, , ]
Fdraw <- Fdraw[thin_val, , ]

#	Average over Gibbs draws
Fdraw2 <- squeeze.mean(Fdraw)
Ldraw2 <- squeeze.mean(Ldraw)
Qdraw2 <- squeeze.mean(Qdraw)
Bdraw2 <- squeeze.mean(Bdraw)

#	Get matrix of autoregressive parameters B
Betas <- NULL
dd <- 1
for(dd in 1:lags){
    beta <- squeeze.mean(Bdraw)
	beta_new <- t(beta)
    if(dd == 1){
		Betas <- beta_new
	  }else{
		Betas <- cbind(Betas, beta_new)
	  }
  }
#
# -- End of "dynfactor_demo.R" --
