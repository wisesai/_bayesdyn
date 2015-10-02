#
#	This program fits a hierarchical dynamic bayesian model averaging. 
#	The output variable is the States' economic growth for Mexico 
#	for the 32 Mexican States. Output is measured using the ITAEE data
#	in forward annual variatios on a quarterly basis.
#
#	Author: 	Juan~Carlos Martínez-Ovando
#	Email:		juan.martinez@banxico.org.mx
#				JC.Martinez.Ovando@gmail.com
#
#	Reference:	Martínez-Ovando, J.~C. (2014) "The Influence of Crime on the Economic Growth 
#					Dynamics in Mexico in the Short-Term," Banco de México, Mimeo.
#
#	Notes:		-	Data must be prepared in Excel before loading it in R 
#					(e.g. growth rates, truncation periods, imputations, etc.)
#				-	Data spans from 2006Q4 to 2011Q3 (i.e. 20 quarterly observations)
#				
#	Updates:	20/12/2014		-	Exporting local data 
#				25/11/2014		-	Revisiting the computation of posterior elasticities and 
#									the computation of posterior model averaging
#				24/09/2014		-	Incorporating informative initial values for state variables
#				15/09/2014		-	Revision of the estimatiom procedure
#				15/09/2014		-	Incorporating new/updated code 
#				13/09/2014		-	Redefining the design matrix
#				28/08/2014		-	Computing different ooutput measurements
#		

rm(list = ls())
rutawork <- 'L:/JCMO.Research/Coding/_r/_my_packages/_bayesdyn/bayesdyn/R/dlmbma'
rutacode <- 'L:/JCMO.Research/Coding/_r/_my_packages/_bayesdyn/bayesdyn/R/dlmbma'
 
#	----------------------------------------------------
#		Loading code
#	----------------------------------------------------
source(paste(rutacode,"/dlmbma_backwardsmoother.R",sep = ""))
source(paste(rutacode,"/dlmbma_demo.R",sep = ""))
source(paste(rutacode,"/dlmbma_forwardfilter.R",sep = ""))
source(paste(rutacode,"/dlmbma_postfilter.R",sep = ""))
source(paste(rutacode,"/dlmbma.R",sep = ""))

#	----------------------------------------------------
#		Loading data
#	----------------------------------------------------
#	Loading the data...
rdata <- read.csv(paste(rutawork,"/xdata.csv",sep = ""), header = TRUE, sep = ",", quote="\"", dec=".", fill = TRUE)

#	Loading matrix of seasonal covariates
Smat <- read.csv(paste(rutawork,"/smat_seasonal.csv",sep = ""), header = TRUE, sep = ",", quote="\"", dec=".", fill = TRUE, row.names = 1)
Smat <- as.matrix(Smat)

Qini <- 20
Qnum <- 20
periodosana <- 	c("2006Q4",
				  "2007Q1","2007Q2","2007Q3","2007Q4",
				  "2008Q1","2008Q2","2008Q3","2008Q4",
				  "2009Q1","2009Q2","2009Q3","2009Q4",
				  "2010Q1","2010Q2","2010Q3","2010Q4",
				  "2011Q1","2011Q2","2011Q3")
  
#	----------------------------------------------------
#		Repository for the output of 32 States
#	----------------------------------------------------

#	Object lists
dma.Xs <- list()

dlmbma.y.filt <- list()
dlmbma.Q.filt <- list()
dlmbma.y.fits <- list()
dlmbma.yQ.fits <- list()
dlmbma.pi.filt <- list()
dlmbma.data.j <- list()
dlmbma.mt.ma <- list()
dlmbma.mt.ma.one <- list()
dlmbma.mt.ma.two <- list()

#	Correlation of adjusted values (MLS estimation, no inference)
rho <- matrix(NA,32,2)
colnames(rho) <- c("Adj.Corr.Y.one","Adj.Corr.Y.two")
rownames(rho) <- c(1:32)

#	----------------------------------------------------
#		Running the model	-	State level
#	----------------------------------------------------
J <- 32 					#	Number of Mexican States
ent <- 1
for(ent in 1:J){
		#'''''''''''''''''''''''''''''''''''''''''''''''		
		#		Setting-up the data set
		#'''''''''''''''''''''''''''''''''''''''''''''''

		#		Indexes
		indexes <- which(rdata[,"clave"] == ent)
		
		#		Annual growth rate (a period ahead)
		Y <- rdata[indexes,c("periodo","itaee_c")]
		Ydat <- as.matrix(Y[Qini:(Qini+Qnum-1),"itaee_c"])
		colnames(Ydat) <- c("itaee_c")
		rownames(Ydat) <- periodosana
		
		#		Covariates X's (initial conditions)
		X <- rdata[indexes,]
		Xmat <- X[Qini:(Qini+Qnum-1), c('itaee_log','itaee_c1','itaee_c2','inpc_general_c1','desocup','ied_mxp_log','usindpro_log','pibp12mes_e','gtopgj_pc')]
		rownames(Xmat) <- periodosana
		
		#		Crime X's 	(criminal incidence)
		Cmat <- X[Qini:(Qini+Qnum-1), c('tclc_log','tflo_log','tfrco_log' )]
		rownames(Cmat) <- periodosana
		
		#		Design matrix	(Note: the function "dyn.dma.R" automatically introduces a constant variable, before S´s, X´s and C´s variables)
		Xdat <- cbind(Smat,Xmat,Cmat)

		Data_j <- cbind(Ydat,Xdat)
		
		Xdat <- as.matrix(Xdat)
		dma.Xs[[ent]] <- Xdat
		Best <- solve(t(Xdat)%*%Xdat)%*%t(Xdat)%*%Ydat
		Yhat.one <- Xdat %*%solve(t(Xdat)%*%Xdat)%*%t(Xdat)%*%Ydat								#	with CRIME variables	
		Yhat.two <- Xdat[,1:13] %*%solve(t(Xdat[,1:13])%*%Xdat[,1:13])%*%t(Xdat[,1:13])%*%Ydat	#	without CRIME variables
		
		rho[ent,1] <- cor(Ydat,Yhat.one) 
		rho[ent,2] <- cor(Ydat,Yhat.two) 

		#	--------------------------------------------
		#		14-09-15	-	Re defining the selection matrix	
		#	--------------------------------------------
		#		Defining which covariables to switch from
		Xwhich <- read.csv(paste(rutawork,"/xmat_which.csv",sep = ""), header = TRUE, sep = ",", quote="\"", dec=".", fill = TRUE)
		modelthreshold <- 21

		#	Executing Dynamic Bayesian Model Averaging
		x <- Xdat[,-1]
		y <- Ydat
		models.which <- Xwhich 
		lambda <- 0.999
		gammaa <- 0.999
		eps <- 0.001/nrow(models.which)

		#'''''''''''''''''''''''''''''''''
		#	Function:	dlmbma
		#'''''''''''''''''''''''''''''''''
		dlmbma_ent <- dlmbma(x, y, models.which, lambda, gammaa, eps, modelthreshold)
		
		#	Saving the output
		#	Filtered
		y.fits <- cbind(Ydat,Yhat.one,Yhat.two,dlmbma_ent$y.filt.ma,dlmbma_ent$y.filt.ma.one,dlmbma_ent$y.filt.ma.two)
		colnames(y.fits) <- c("itaee_c","itaee_yhat.one","itaee_yhat.two","itaee_yfilt_ma","itaee_yfilt_ma.one","itaee_yfilt_ma.two")

		#	Predictive dispersion
		yQ.fits <- cbind(Ydat,dlmbma_ent$y.filt.ma,sqrt(dlmbma_ent$Q.filt.ma),dlmbma_ent$y.filt.ma.one,sqrt(dlmbma_ent$Q.filt.ma.one),dlmbma_ent$y.filt.ma.two,sqrt(dlmbma_ent$Q.filt.ma.two))
		colnames(yQ.fits) <- c("itaee_c","itaee_yfilt_ma","itaee_Qfilt_ma","itaee_yfilt_ma.one","itaee_Qfilt_ma.one","itaee_yfilt_ma.two","itaee_Qfilt_ma.two")
		
		dlmbma.y.filt[[ent]] <- dlmbma_ent$y.filt
		dlmbma.Q.filt[[ent]] <- dlmbma_ent$Q.filt
		dlmbma.y.fits[[ent]] <- y.fits
		dlmbma.yQ.fits[[ent]] <- yQ.fits
		dlmbma.pi.filt[[ent]] <- dlmbma_ent$pi.filt
		dlmbma.data.j[[ent]] <- Data_j
		dlmbma.mt.ma[[ent]] <- dlmbma_ent$m.filt.ma
		dlmbma.mt.ma.one[[ent]] <- dlmbma_ent$m.filt.ma.one
		dlmbma.mt.ma.two[[ent]] <- dlmbma_ent$m.filt.ma.two
		
	# -- END of "dlmbma.R" for the 32 States --
	}
	
#''''''''''''''''''''''''''''''''''''''''''''''''
#	Exporting CSV outputs
#''''''''''''''''''''''''''''''''''''''''''''''''
#	Exporting the mt's for general BMA and scenarios ONE and TWO
ent <- 1
for(ent in 1:J){
	matrix.aux <- dlmbma.mt.ma[[ent]]
	write.table(matrix.aux, 
				file = paste(rutawork,"/output_dlmbma/mt_",ent,".csv",sep = ""), sep = ",", col.names = TRUE, qmethod = "double")
	matrix.aux <- dlmbma.mt.ma.one[[ent]]
	write.table(matrix.aux, 
				file = paste(rutawork,"/output_dlmbma/mt_",ent,"_one.csv",sep = ""), sep = ",", col.names = TRUE, qmethod = "double")
	matrix.aux <- dlmbma.mt.ma.two[[ent]]
	write.table(matrix.aux, 
				file = paste(rutawork,"/output_dlmbma/mt_",ent,"_two.csv",sep = ""), sep = ",", col.names = TRUE, qmethod = "double")
	}
#	Exporting Data_j
ent <- 1
for(ent in 1:J){
	matrix.aux <- dlmbma.data.j[[ent]]
	write.table(matrix.aux, 
				file = paste(rutawork,"/output_dlmbma/data_",ent,".csv",sep = ""), sep = ",", col.names = TRUE, qmethod = "double")
	}
#	Exporting posterior filtering paths
ent <- 1
for(ent in 1:J){
	matrix.aux <- dlmbma.y.filt[[ent]]
	write.table(matrix.aux, 
				file = paste(rutawork,"/output_dlmbma/filter_",ent,".csv",sep = ""), sep = ",", col.names = TRUE, qmethod = "double")
	matrix.aux <- dlmbma.Q.filt[[ent]]
	write.table(matrix.aux, 
				file = paste(rutawork,"/output_dlmbma/filter_",ent,"_Q.csv",sep = ""), sep = ",", col.names = TRUE, qmethod = "double")
	}
#	Exporting posterior filtering paths
ent <- 1
for(ent in 1:J){
	matrix.aux <- dlmbma.y.fits[[ent]]
	write.table(matrix.aux, 
				file = paste(rutawork,"/output_dlmbma/fit_",ent,".csv",sep = ""), sep = ",", col.names = TRUE, qmethod = "double")
	matrix.aux <- dlmbma.yQ.fits[[ent]]
	write.table(matrix.aux, 
				file = paste(rutawork,"/output_dlmbma/fit_",ent,"_Q.csv",sep = ""), sep = ",", col.names = TRUE, qmethod = "double")
	}
#	Exporting posterior probabilities for the alternative models
ent <- 1
for(ent in 1:J){
	matrix.aux <- dlmbma.pi.filt[[ent]]
	write.table(matrix.aux, 
				file = paste(rutawork,"/output_dlmbma/pmp_",ent,".csv",sep = ""), sep = ",", col.names = TRUE, qmethod = "double")
	}
	
#	----------------------------------------------------
#		Loading additional data
#	----------------------------------------------------
#	Loading the weights to recover the aggregated ITAEN
pib_pesos <- read.csv(paste(rutawork,"/smat_weitghs.csv",sep = ""), header = TRUE, sep = ",", quote="\"", dec=".", fill = TRUE)
pib_pesos <- as.data.frame(pib_pesos)

#	Loading matching variables for States and Regions (Banco de México)
claves_regiones <- read.csv(paste(rutawork,"/sclaves_regiones.csv",sep = ""), header = TRUE, sep = ",", quote="\"", dec=".", fill = TRUE)
claves_regiones <- as.data.frame(claves_regiones)

#========================================================================
#		Inference	-	Regional level
#========================================================================

	#************************************************************************
	#			Yfilt.ma	-	Regiones
	#************************************************************************

	TT <- 20

	#----------------------------------
	# 	Region NORTE
	#----------------------------------
	ents_region <- claves_regiones[which(claves_regiones[,"clavereg"]==1),"clave"] 
	N_region <- length(ents_region)
	y.filt.region <- array(NA,dim=c(TT,27,N_region))
	pi.filt.region <- array(NA,dim=c(TT,27,N_region))
	y.filt.models <- array(NA,dim=c(TT,27))
	colnames(y.filt.models) <- c(1:27)
	rownames(y.filt.models) <- periodosana
	y.act.data <- array(NA,dim=c(TT,N_region))
	colnames(y.act.data) <- ents_region
	rownames(y.act.data) <- periodosana
	y.fit.data <- array(NA,dim=c(TT,N_region))
	colnames(y.fit.data) <- ents_region
	rownames(y.fit.data) <- periodosana
	pi.filt.crime <- array(NA,dim=c(TT,3))
	colnames(pi.filt.crime) <- c('crime','no.crime','post_odd')
	rownames(pi.filt.crime) <- periodosana
	#	States' shares for the region
	pib_region <- t(as.matrix(pib_pesos[ents_region,"pib"])) / sum(as.matrix(pib_pesos[ents_region,"pib"]))
	rownames(pib_region) <- c("pib.norte")
	colnames(pib_region) <- c(ents_region)
	#	Nested model averaging
	n <- 1 
	for(n in 1:N_region){
		y.filt.region[,,n] <- dlmbma.y.filt[[ents_region[n]]]
		pi.filt.region[,,n] <- dlmbma.pi.filt[[ents_region[n]]]
		}
	#	Aggregate nested model averaging
	j <- 1
	for(j in 1:27){
		y.filt.models[,j] <- (y.filt.region[,j,] * pi.filt.region[,j,]) %*% t(pib_region)
		}
	y.filt.norte <- as.matrix(rowSums(y.filt.models))
	colnames(y.filt.norte) <- c("itaee_filt_norte")
	rownames(y.filt.norte) <- periodosana
	#	Aggregate actual and fitted data
	ent <- 1
	for(ent in 1:N_region){	
		#		Annual growth rate (a period ahead)
		y.act.data[,ent] <- dlmbma.y.fits[[ents_region[ent]]][,"itaee_c"]
		y.fit.data[,ent] <- dlmbma.y.fits[[ents_region[ent]]][,"itaee_yhat.one"]
		}
	y.dat.norte <- y.act.data %*% t(pib_region)	
	colnames(y.dat.norte) <- c("itaee_c_norte")
	y.fit.norte <- y.fit.data %*% t(pib_region)	
	colnames(y.fit.norte) <- c("itaee_c_norte")
	y.fit.r.norte <- cbind(y.dat.norte,y.filt.norte,y.fit.norte)
	write.table(y.fit.r.norte, 
				file = paste(rutawork,"/output_dlmbma/fit_region_norte.csv",sep = ""), sep = ",", col.names = TRUE, qmethod = "double")

	#----------------------------------
	# 	Region CENTRO-NORTE
	#----------------------------------
	ents_region <- claves_regiones[which(claves_regiones[,"clavereg"]==2),"clave"] 
	N_region <- length(ents_region)
	y.filt.region <- array(NA,dim=c(TT,27,N_region))
	pi.filt.region <- array(NA,dim=c(TT,27,N_region))
	y.filt.models <- array(NA,dim=c(TT,27))
	colnames(y.filt.models) <- c(1:27)
	rownames(y.filt.models) <- periodosana
	y.act.data <- array(NA,dim=c(TT,N_region))
	colnames(y.act.data) <- ents_region
	rownames(y.act.data) <- periodosana
	y.fit.data <- array(NA,dim=c(TT,N_region))
	colnames(y.fit.data) <- ents_region
	rownames(y.fit.data) <- periodosana
	#	States' shares for the region
	pib_region <- t(as.matrix(pib_pesos[ents_region,"pib"])) / sum(as.matrix(pib_pesos[ents_region,"pib"]))
	rownames(pib_region) <- c("pib.centronorte")
	colnames(pib_region) <- c(ents_region)
	#	Nested model averaging
	n <- 1 
	for(n in 1:N_region){
		y.filt.region[,,n] <- dlmbma.y.filt[[ents_region[n]]]
		pi.filt.region[,,n] <- dlmbma.pi.filt[[ents_region[n]]]
		}
	#	Aggregate nested model averaging
	j <- 1
	for(j in 1:27){
		y.filt.models[,j] <- (y.filt.region[,j,] * pi.filt.region[,j,]) %*% t(pib_region)
		}
	y.filt.centronorte <- as.matrix(rowSums(y.filt.models))
	colnames(y.filt.centronorte) <- c("itaee_filt_centronorte")
	rownames(y.filt.centronorte) <- periodosana
	#	Aggregate actual and fitted data
	ent <- 1
	for(ent in 1:N_region){	
		#		Annual growth rate (a period ahead)
		y.act.data[,ent] <- dlmbma.y.fits[[ents_region[ent]]][,"itaee_c"]
		y.fit.data[,ent] <- dlmbma.y.fits[[ents_region[ent]]][,"itaee_yhat.one"]
		}
	y.dat.centronorte <- y.act.data %*% t(pib_region)	
	colnames(y.dat.centronorte) <- c("itaee_c_centronorte")
	y.fit.centronorte <- y.fit.data %*% t(pib_region)	
	colnames(y.fit.centronorte) <- c("itaee_fit_centronorte")
	y.fit.r.centronorte <- cbind(y.dat.centronorte,y.filt.centronorte,y.fit.centronorte)
	write.table(y.fit.r.centronorte, 
				file = paste(rutawork,"/output_dlmbma/fit_region_centronorte.csv",sep = ""), sep = ",", col.names = TRUE, qmethod = "double")
	
	#----------------------------------
	# 	Region CENTRO
	#----------------------------------
	ents_region <- claves_regiones[which(claves_regiones[,"clavereg"]==3),"clave"] 
	N_region <- length(ents_region)
	y.filt.region <- array(NA,dim=c(TT,27,N_region))
	pi.filt.region <- array(NA,dim=c(TT,27,N_region))
	y.filt.models <- array(NA,dim=c(TT,27))
	colnames(y.filt.models) <- c(1:27)
	rownames(y.filt.models) <- periodosana
	y.act.data <- array(NA,dim=c(TT,N_region))
	colnames(y.act.data) <- ents_region
	rownames(y.act.data) <- periodosana
	y.fit.data <- array(NA,dim=c(TT,N_region))
	colnames(y.fit.data) <- ents_region
	rownames(y.fit.data) <- periodosana
	#	States' shares for the region
	pib_region <- t(as.matrix(pib_pesos[ents_region,"pib"])) / sum(as.matrix(pib_pesos[ents_region,"pib"]))
	rownames(pib_region) <- c("pib.centro")
	colnames(pib_region) <- c(ents_region)
	#	Nested model averaging
	n <- 1 
	for(n in 1:N_region){
		y.filt.region[,,n] <- dlmbma.y.filt[[ents_region[n]]]
		pi.filt.region[,,n] <- dlmbma.pi.filt[[ents_region[n]]]
		}
	#	Aggregate nested model averaging
	j <- 1
	for(j in 1:27){
		y.filt.models[,j] <- (y.filt.region[,j,] * pi.filt.region[,j,]) %*% t(pib_region)
		}
	y.filt.centro <- as.matrix(rowSums(y.filt.models))
	colnames(y.filt.centro) <- c("itaee_filt_centro")
	rownames(y.filt.centro) <- periodosana
	#	Aggregate actual and fitted data
	ent <- 1
	for(ent in 1:N_region){	
		#		Annual growth rate (a period ahead)
		y.act.data[,ent] <- dlmbma.y.fits[[ents_region[ent]]][,"itaee_c"]
		y.fit.data[,ent] <- dlmbma.y.fits[[ents_region[ent]]][,"itaee_yhat.one"]
		}
	y.dat.centro <- y.act.data %*% t(pib_region)	
	colnames(y.dat.centro) <- c("itaee_c_centro")
	y.fit.centro <- y.fit.data %*% t(pib_region)	
	colnames(y.fit.centro) <- c("itaee_fit_centro")
	y.fit.r.centro <- cbind(y.dat.centro,y.filt.centro,y.fit.centro)
	write.table(y.fit.r.centro, 
				file = paste(rutawork,"/output_dlmbma/fit_region_centro.csv",sep = ""), sep = ",", col.names = TRUE, qmethod = "double")
	
	#----------------------------------
	# 	Region SUR
	#----------------------------------
	ents_region <- claves_regiones[which(claves_regiones[,"clavereg"]==4),"clave"] 
	N_region <- length(ents_region)
	y.filt.region <- array(NA,dim=c(TT,27,N_region))
	pi.filt.region <- array(NA,dim=c(TT,27,N_region))
	y.filt.models <- array(NA,dim=c(TT,27))
	colnames(y.filt.models) <- c(1:27)
	rownames(y.filt.models) <- periodosana
	y.act.data <- array(NA,dim=c(TT,N_region))
	colnames(y.act.data) <- ents_region
	rownames(y.act.data) <- periodosana
	y.fit.data <- array(NA,dim=c(TT,N_region))
	colnames(y.fit.data) <- ents_region
	rownames(y.fit.data) <- periodosana
	#	States' shares for the region
	pib_region <- t(as.matrix(pib_pesos[ents_region,"pib"])) / sum(as.matrix(pib_pesos[ents_region,"pib"]))
	rownames(pib_region) <- c("pib.sur")
	colnames(pib_region) <- c(ents_region)
	#	Nested model averaging
	n <- 1 
	for(n in 1:N_region){
		y.filt.region[,,n] <- dlmbma.y.filt[[ents_region[n]]]
		pi.filt.region[,,n] <- dlmbma.pi.filt[[ents_region[n]]]
		}
	#	Aggregate nested model averaging
	j <- 1
	for(j in 1:27){
		y.filt.models[,j] <- (y.filt.region[,j,] * pi.filt.region[,j,]) %*% t(pib_region)
		}
	y.filt.sur <- as.matrix(rowSums(y.filt.models))
	colnames(y.filt.sur) <- c("itaee_filt_sur")
	rownames(y.filt.sur) <- periodosana
	#	Aggregate actual and fitted data
	ent <- 1
	for(ent in 1:N_region){	
		#		Annual growth rate (a period ahead)
		y.act.data[,ent] <- dlmbma.y.fits[[ents_region[ent]]][,"itaee_c"]
		y.fit.data[,ent] <- dlmbma.y.fits[[ents_region[ent]]][,"itaee_yhat.one"]
		}
	y.dat.sur <- y.act.data %*% t(pib_region)	
	colnames(y.dat.sur) <- c("itaee_c_sur")
	y.fit.sur <- y.fit.data %*% t(pib_region)	
	colnames(y.fit.sur) <- c("itaee_fit_sur")
	y.fit.r.sur <- cbind(y.dat.sur,y.filt.sur,y.fit.sur)
	write.table(y.fit.r.sur, 
				file = paste(rutawork,"/output_dlmbma/fit_region_sur.csv",sep = ""), sep = ",", col.names = TRUE, qmethod = "double")


	#----------------------------------
	# 	Region NATIONAL
	#----------------------------------
	ents_region <- claves_regiones[ which( (claves_regiones[,"clavereg"]==1) | 
										   (claves_regiones[,"clavereg"]==2) | 
										   (claves_regiones[,"clavereg"]==3) | 
										   (claves_regiones[,"clavereg"]==4) ) ,"clave"] 
	N_region <- length(ents_region)
	y.filt.region <- array(NA,dim=c(TT,27,N_region))
	pi.filt.region <- array(NA,dim=c(TT,27,N_region))
	y.filt.models <- array(NA,dim=c(TT,27))
	colnames(y.filt.models) <- c(1:27)
	rownames(y.filt.models) <- periodosana
	y.act.data <- array(NA,dim=c(TT,N_region))
	colnames(y.act.data) <- ents_region
	rownames(y.act.data) <- periodosana
	y.fit.data <- array(NA,dim=c(TT,N_region))
	colnames(y.fit.data) <- ents_region
	rownames(y.fit.data) <- periodosana
	#	States' shares for the region
	pib_region <- t(as.matrix(pib_pesos[ents_region,"pib"])) / sum(as.matrix(pib_pesos[ents_region,"pib"]))
	rownames(pib_region) <- c("pib.nacional")
	colnames(pib_region) <- c(ents_region)
	#	Nested model averaging
	n <- 1 
	for(n in 1:N_region){
		y.filt.region[,,n] <- dlmbma.y.filt[[ents_region[n]]]
		pi.filt.region[,,n] <- dlmbma.pi.filt[[ents_region[n]]]
		}
	#	Aggregate nested model averaging
	j <- 1
	for(j in 1:27){
		y.filt.models[,j] <- (y.filt.region[,j,] * pi.filt.region[,j,]) %*% t(pib_region)
		}
	y.filt.nacional <- as.matrix(rowSums(y.filt.models))
	colnames(y.filt.nacional) <- c("itaee_filt_nacional")
	rownames(y.filt.nacional) <- periodosana
	#	Aggregate actual and fitted data
	ent <- 1
	for(ent in 1:N_region){	
		#		Annual growth rate (a period ahead)
		y.act.data[,ent] <- dlmbma.y.fits[[ents_region[ent]]][,"itaee_c"]
		y.fit.data[,ent] <- dlmbma.y.fits[[ents_region[ent]]][,"itaee_yhat.one"]
		}
	y.dat.nacional <- y.act.data %*% t(pib_region)	
	colnames(y.dat.nacional) <- c("itaee_c_nacional")
	y.fit.nacional <- y.fit.data %*% t(pib_region)	
	colnames(y.fit.nacional) <- c("itaee_fit_nacional")
	y.fit.r.nacional <- cbind(y.dat.nacional,y.filt.nacional,y.fit.nacional)
	write.table(y.fit.r.nacional, 
				file = paste(rutawork,"/output_dlmbma/fit_region_nacional.csv",sep = ""), sep = ",", col.names = TRUE, qmethod = "double")


#	Saving outputs
save( dlmbma.y.filt, dlmbma.Q.filt,
	  dlmbma.y.fits, dlmbma.yQ.fits,
	  dlmbma.pi.filt,
	  dlmbma.data.j,
	  dlmbma.mt.ma,
	  dlmbma.mt.ma.one,
	  dlmbma.mt.ma.two,
	  file = paste(rutawork,"/output_dlmbma/dlmbma_out.RData",sep = "")
	)
	
#
#		END of "dlmbma_demo.R"	