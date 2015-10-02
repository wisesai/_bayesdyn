#
#   tvpregsv.R - demo
#
#   Author: 	Juan Carlos Martinez-Ovando 
#				(Email: juan.martinez@banxico.org.mx)
#				(Email: JC.Martinez.Ovando@gmail.com)
#
#	Revision:		May 17, 2015
#

rm(list = ls())

#	Packages
#install.packages('bnormnlr')
library('bnormnlr')
#install.packages('MASS')
library('MASS')

#	Paths
path.code <- 'L:/JCMO.Research/Banxico/proyectos/2014/_ctot&inflationmx/_code_r'
path.data <- 'L:/JCMO.Research/Banxico/proyectos/2014/_ctot&inflationmx/_datos/_datos_tesis_juan'
path.output <- 'L:/JCMO.Research/Banxico/proyectos/2014/_ctot&inflationmx/_code_r'

#	Source code
source(paste(path.code,"/rowprodmatrix.R",sep=""))
source(paste(path.code,"/uni.slice.R",sep=""))
source(paste(path.code,"/tvpregsv_beta.R",sep=""))
source(paste(path.code,"/tvpregsv_alpha.R",sep=""))
source(paste(path.code,"/tvpregsv_gamma.R",sep=""))
source(paste(path.code,"/tvpregsv_h.R",sep=""))
source(paste(path.code,"/tvpregsv_eta.R",sep=""))
source(paste(path.code,"/tvpregsv_varphi.R",sep=""))
source(paste(path.code,"/tvpregsv_phi.R",sep=""))
source(paste(path.code,"/tvpregsv_sigmah.R",sep=""))
source(paste(path.code,"/tvpregsv.R",sep=""))

#	Data
inpc <- read.csv( paste(path.data,"/banxico_inpc_no_subyacente.csv",sep=""), row.names = 1, header = TRUE )
ctot <- read.csv( paste(path.data,"/Commodity_Terms_of_Trade_MX_150312.csv",sep=""), row.names = 1, header = TRUE )

y <- as.matrix(inpc[,"inf_ns_a_fv"])
x <- cbind( matrix(1,dim(y)[1],1), as.matrix(inpc[,"inf_ns"]))
w <- as.matrix(ctot[rownames(inpc),"ctot_mx"])
ctot[ ,"vctot_mx"] <- (as.matrix(ctot[,"ctot_mx"])-mean(ctot[,"ctot_mx"]))^2
z <- as.matrix(ctot[rownames(inpc),"vctot_mx"])

N_mcmc <- 3
N_slice <- 13

uni.slice.calls <<- 0	# Number of calls of the slice sampling function
uni.slice.evals <<- 0	# Number of density evaluations done in these calls

#
#	MCMC sampling algorithm
#
tvpregsv_out <- tvpregsv(y,x,w,z,N_mcmc,N_slice)
beta_mcmc <- tvpregsv_out[[1]]
alpha_mcmc <- tvpregsv_out[[2]]
gamma_mcmc <- tvpregsv_out[[3]]
h_mcmc <- tvpregsv_out[[4]]
eta_mcmc <- tvpregsv_out[[5]]
varphi_mcmc <- tvpregsv_out[[6]]
phi_mcmc <- tvpregsv_out[[7]]
sigmah_mcmc <- tvpregsv_out[[8]]

save(tvpregsv_out, file = paste(path.output,"/tvpregsv_out.Rdata",sep=""))
 
#
#	--	END	--