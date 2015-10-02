#
#	Demo for "tvpvar.R"
#

rm(list=ls())

#	Set working directory
#	setwd("D:/Proyectos/R_DB")
ruta.dir <- "L:/JCMO.Research/Coding/_r/_my_packages/_bayesdyn/bayesdyn/R/tvpvar"
ruta.code <- "L:/JCMO.Research/Coding/_r/_my_packages/_bayesdyn/bayesdyn/R"
setwd(ruta.dir)

#	Packages
#install.packages('mvtnorm')
library('mvtnorm')

#	Source code
source(paste(ruta.code,"/utilities/colStdev.R",sep=""))
source(paste(ruta.code,"/utilities/corrvc.R",sep=""))
source(paste(ruta.code,"/utilities/mlag2.R",sep=""))
source(paste(ruta.code,"/utilities/reshaper.R",sep=""))
source(paste(ruta.code,"/utilities/reshapes.R",sep=""))
source(paste(ruta.code,"/utilities/squeeze.sdev.R",sep=""))
source(paste(ruta.code,"/utilities/wishartrnd.R",sep=""))

source(paste(ruta.code,"/tvpvar/tvpvar_prior.R",sep=""))
source(paste(ruta.code,"/tvpvar/tvpvar_primiceri.R",sep=""))
source(paste(ruta.code,"/tvpvar/gibbs_ss.R",sep=""))
source(paste(ruta.code,"/tvpvar/gibbs_B.R",sep=""))
source(paste(ruta.code,"/tvpvar/gibbs_ssa.R",sep=""))
source(paste(ruta.code,"/tvpvar/gibbs_Q.R",sep=""))
source(paste(ruta.code,"/tvpvar/gibbs_A.R",sep=""))
source(paste(ruta.code,"/tvpvar/gibbs_logSigma.R",sep=""))
source(paste(ruta.code,"/tvpvar/gibbs_ssb.R",sep=""))
source(paste(ruta.code,"/tvpvar/tvpvar.R",sep=""))

#	Reading data file from Matlab
ydata <- read.csv(paste(ruta.dir,"/ydata.csv",sep=""), header = TRUE, row.names = 1)
ydata <- as.matrix(ydata)

#	Number of factors & lags:
N_training_prior <- 40		#	N_training_prior is the size of the training sample
p_lags <- 2					# 	p_lags is number of lags in the VAR part

# 	Set some Gibbs - related preliminaries
N_mcmc <- 100				#	Number of replications
N_burn <- 20				#	Number of burn-in-draws
N_print <- 10				#	Print in the screen every "N_print"-th iteration

prior_type <- 1				#	informative prior = 1, non informative prior = 0
demean <- 0					#	standardize the data = 1; otherwise = 0

#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#									Gibbs sampler
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
tvpvar_out <- tvpvar(ydata, demean, N_training_prior, p_lags, N_mcmc, N_burn, N_print, prior_type)
Bt_postmean <- tvpvar_out[[1]]
At_postmean <- tvpvar_out[[2]]
Sigt_postmean <- tvpvar_out[[3]]
Qmean <- tvpvar_out[[4]]
Smean <- tvpvar_out[[5]]
Wmean <- tvpvar_out[[6]]
sigmean <- tvpvar_out[[7]]
cormean <- tvpvar_out[[8]]

Bt_postmean <- Bt_postmean / N_mcmc								#	Posterior mean of B(t) (VAR regression coeff.)
At_postmean <- At_postmean / N_mcmc								#	Posterior mean of A(t) (VAR covariances)
Sigt_postmean <- Sigt_postmean / N_mcmc							#	Posterior mean of SIGMA(t) (VAR variances)
Qmean <- Qmean / N_mcmc											#	Posterior mean of Q (covariance of B(t))
Smean <- Smean / N_mcmc											#	Posterior mean of S (covariance of A(t))
Wmean <- Wmean / N_mcmc											#	Posterior mean of W (covariance of SIGMA(t))

sigmean <- sigmean / N_mcmc
cormean <- cormean / N_mcmc 
sig2mo <- sig2mo / N_mcmc
cor2mo <- cor2mo / N_mcmc

#	Standard deviations of residuals of Inflation, Unemployment and Interest Rates
pdf(paste(ruta.dir,"/plot_sdresiduals_primiceri.pdf",sep=""))
par(mfrow=c(2,2))
colnames(sigmean) <- c("Inflation", "Unemplyment", "Interest Rate")
sigmean_ts <- ts(sigmean, frequency = 4, start = c(1963, 2))
plot(sigmean_ts, xlab = "Quarter", ylab = "Residual", main = "Posterior mean of the standard deviation of residuals") 
dev.off()

#	Output
tvpvar_data <- list(ydata, yearlab)
tvpvar_mcmc <- list(N_mcmc, N_burn, N_print)
tvpvar_prior <- list(B_OLS, VB_OLS, A_OLS, sigma_OLS, VA_OLS)
tvpvar_output <- list(Bt_postmean, At_postmean, Sigt_postmean, Qmean, Smean, Wmean, sigmean, cormean, sig2mo, cor2mo)

#	Saving output...
save(tvpvar_data, tvpvar_mcmc, tvpvar_prior, tvpvar_output, file = paste(ruta.dir, "/tvpvar_output_primiceri.RData", sep=""))

#
#	-- End of "Hetero_TVP_VAR.R" --