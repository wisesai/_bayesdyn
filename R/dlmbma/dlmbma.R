dlmbma <- function(x, y, models.which, lambda=0.99, gammaa=0.99, eps=.001/nrow(models.which), modelthreshold=21){
#
#	dlmbma.R
#	Bayesian Dynamic Model Averaging (Filtering and Smoothing)
#	This function computes Bayesian model avraging, forward-filtering and backward smoothing of a general 
#	homoscedastic regression dynamic linear model among several model confiburations.
#
# 	INPUTS:
# 		x					-	TTxp matrix of system inputs
# 		y					-	TT-vector of system outputs
# 		models.which		-	Kxp matrix, with 1 row per model and 1 col per variable	indicating whether that 
#								variable is in the model (the state theta is of dim (model.dim+1); the extra 1 for the intercept)
# 		lambda				-	Discount factor (W&H's notion)
# 		gammaa				-	Flatterning parameter for model updating
# 		eps					-	Regularization parameter for regularizing posterior model probabilities away from zero
#		modelthreshold 		-	Threshold for the model including and excluing "causal" variables.
#
# 	OUTPUTS:
# 		yhat.bymodel		-	TTxK matrix whose tk element gives yhat for yt for model k
# 		yhat.ma				-	TT vector whose t element gives the model-averaged yhat for yt
# 		yhat.ma.one			-	TT vector whose t element gives the model-averaged yhat for yt for models {1,..., "modelthreshold"}
# 		yhat.ma.two			-	TT vector whose t element gives the model-averaged yhat for yt for models {"modelthreshold"+1, ..., K}
# 		pmp					-	TTxK matrix whose tk element is the post prob of model k at t
# 		thetahat.ma			-	TTx(nvar+1) matrix whose tk element is the model-averaged estimate of theta_{j-1} at t
# 		Vtheta.ma			-	TTx(nvar+1) matrix whose tk element is the model-average variance of thetahat_{j-1} at t
#
#	Author: 	Juan~Carlos Martínez-Ovando
#	Email:		juan.martinez@banxico.org.mx
#				jc.martinez.ovando@gmail.com
#
#	Reference:	Martínez-Ovando, J.~C. (2014) "The Influence of Crime on the Economic Growth 
#					Dynamics in Mexico in the Short-Term," Banco de México, Mimeo.
#

#	Initialization
y <- as.matrix(y)
x <- as.matrix(x)
nvar <- ncol(x)									# 	Number of variables	
TT <- length(y)									#	Time span
KK <- nrow (models.which)						#	Number of models
if( (modelthreshold < 2) | (modelthreshold > (KK-1)) ){print("Error in the specification of the model threshold for the two scenarios."); break}
model.dim <- rep(0,KK)							#	Dimensions of the models (according to Xwhich)
for(k in 1:KK){
	model.dim[k] <- sum(as.real(models.which[k,]))+1
	}
theta.filt <- vector("list",KK)					#	Object list for filtered $\theta$
Sigma.filt <- vector("list",KK) 				#	Object list for filtered variance for $\theta$

m_models <- matrix(0,(nvar+1),KK)				#	Repository for the filtered 'm_t's
rownames(m_models) <- c('cte',colnames(x))
colnames(m_models) <- c(1:KK)

m.filt.ma <- matrix(0,TT,(nvar+1))				#	Repository for m_t´s BMA
rownames(m.filt.ma) <- rownames(x) 
colnames(m.filt.ma) <- c('cte',colnames(x))
m.filt.ma.one <- m.filt.ma						#	Repository for m_t´s BMA scenario ONE
m.filt.ma.two <- m.filt.ma						#	Repository for m_t´s BMA scenario TWO

theta.smooth <- vector ("list",KK)				#	Object list for smoothed $\theta$
Sigma.smooth <- vector ("list",KK) 				#	Object list for smoothed variance for $\theta$

a.filt <- vector("list",TT*KK)					#	Object list for filtered $a$
R.filt <- vector("list",TT*KK) 					#	Object list for filtered variance for $R$

y.filt <- matrix (rep(0,TT*KK), ncol=KK)		#	Repository for Y filtered
y.smooth <- matrix (rep(0,TT*KK), ncol=KK)		#	Repository for Y smoothed
rownames(y.filt) <- rownames(y)
rownames(y.smooth) <- rownames(y)
colnames(y.filt) <- c(1:KK)
colnames(y.smooth) <- c(1:KK)

Q.filt <- matrix (rep(0,TT*KK), ncol=KK)		#	Repository for the variance of Y filtered
Q.smooth <- matrix (rep(0,TT*KK), ncol=KK)		#	Repository for the variance of Y smoothed
rownames(Q.filt) <- rownames(y)
rownames(Q.smooth) <- rownames(y)
colnames(Q.filt) <- c(1:KK)
colnames(Q.smooth) <- c(1:KK)

pi.filt <- matrix (rep(0,TT*KK), ncol=KK)		#	Repository for posterior model probabilities filtered
pi.smooth <- matrix (rep(0,TT*KK), ncol=KK)		#	Repository for posterior model probabilities smoothed
rownames(pi.filt) <- rownames(y)
rownames(pi.smooth) <- rownames(y)
colnames(pi.filt) <- c(1:KK)
colnames(pi.smooth) <- c(1:KK)

y.filt.ma <- as.matrix(rep(0,TT))				#	Repository for Y filtered BMA
y.smooth.ma <- as.matrix(rep(0,TT))				#	Repository for Y smoothed BMA
rownames(y.filt.ma) <- rownames(y)
rownames(y.smooth.ma) <- rownames(y)
colnames(y.filt.ma) <- c("y.filt.ma")
colnames(y.smooth.ma) <- c("y.smooth.ma")

Q.filt.ma <- as.matrix(rep(0,TT))				#	Repository for the variance Y filtered BMA
Q.smooth.ma <- as.matrix(rep(0,TT))				#	Repository for the variance Y smoothed BMA
rownames(Q.filt.ma) <- rownames(y)
rownames(Q.smooth.ma) <- rownames(y)
colnames(Q.filt.ma) <- c("Q.filt.ma")
colnames(Q.smooth.ma) <- c("Q.smooth.ma")

y.filt.ma.one <- as.matrix(rep(0,TT))			#	Repository for Y filtered BMA scenario I
y.smooth.ma.one <- as.matrix(rep(0,TT))			#	Repository for Y smoothed BMA scenario I
rownames(y.filt.ma.one) <- rownames(y)
rownames(y.smooth.ma.one) <- rownames(y)
colnames(y.filt.ma.one) <- c("y.filt.ma.one")
colnames(y.smooth.ma.one) <- c("y.smooth.ma.one")

Q.filt.ma.one <- as.matrix(rep(0,TT))			#	Repository for the variance Y filtered BMA scenario I
Q.smooth.ma.one <- as.matrix(rep(0,TT))			#	Repository for the variance Y smoothed BMA scenario I
rownames(Q.filt.ma.one) <- rownames(y)
rownames(Q.smooth.ma.one) <- rownames(y)
colnames(Q.filt.ma.one) <- c("Q.filt.ma.one")
colnames(Q.smooth.ma.one) <- c("Q.smooth.ma.one")

y.filt.ma.two <- as.matrix(rep(0,TT))			#	Repository for Y filtered BMA scenario II
y.smooth.ma.two <- as.matrix(rep(0,TT))			#	Repository for Y smoothed BMA scenario II
rownames(y.filt.ma.two) <- rownames(y)
rownames(y.smooth.ma.two) <- rownames(y)
colnames(y.filt.ma.two) <- c("y.filt.ma.two")
colnames(y.smooth.ma.two) <- c("y.smooth.ma.two")

Q.filt.ma.two <- as.matrix(rep(0,TT))			#	Repository for the variance Y filtered BMA scenario II
Q.smooth.ma.two <- as.matrix(rep(0,TT))			#	Repository for the variance Y smoothed BMA scenario II
rownames(Q.filt.ma.two) <- rownames(y)
rownames(Q.smooth.ma.two) <- rownames(y)
colnames(Q.filt.ma.two) <- c("Q.filt.ma.two")
colnames(Q.smooth.ma.two) <- c("Q.smooth.ma.two")

#'''''''''''''''''''''''''''''''''''''''''''
#	Additional variables
#'''''''''''''''''''''''''''''''''''''''''''
cte <- as.matrix(rep(1,TT))						#	Constant term
colnames(cte) <- c("cte")
rownames(cte) <- rownames(y)

#	Temporaty repositories
m0.k <- vector("list",KK)
C0.k <- vector("list",KK)
n0.k <- vector("list",KK)
s0.k <- vector("list",KK)

#'''''''''''''''''''''''''''''''''''''''''''
# 	Iterations
#'''''''''''''''''''''''''''''''''''''''''''
#	Model Threshold (Two Scenarios Only)
modelthreshold <- 21
#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#	Forward Filtering
#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
tt <- 1	
for(tt in 1:TT){	
	kk <- 1
	for (kk in 1:KK){
		#	Reading covariates
		which.Xs <- as.logical( models.which[kk,])
		#	Automatically includes the constant term
		Xwhich.k <- cbind(cte,x[,which.Xs])
		nvar.k <- ncol(Xwhich.k)
		#
		#	Initializing filtering variables
		#
		if(tt == 1){
			m0.k[[kk]] <- solve(t(Xwhich.k)%*%Xwhich.k)%*%t(Xwhich.k)%*%y
			e0 <- y - Xwhich.k%*%solve(t(Xwhich.k)%*%Xwhich.k)%*%t(Xwhich.k)%*%y 
			C0.k[[kk]] <- solve(t(Xwhich.k)%*%Xwhich.k)
			n0.k[[kk]] <- TT
			s0.k[[kk]] <- (t(e0)%*%e0)/(n0.k[[kk]]-1)	
			}else{
			}
		m0 <- m0.k[[kk]]
		C0 <- C0.k[[kk]]
		n0 <- n0.k[[kk]]
		s0 <- s0.k[[kk]]
		Ft <- Xwhich.k[tt,]
		Gt <- diag(nvar.k)
		yt <- y[tt]
		mdl.filtered <- dlmbma_forwardfilter(m0,C0,n0,s0,Ft,Gt,yt)
		
		#	Updating initial conditions
		m0.k[[kk]] <- mdl.filtered$mt
		C0.k[[kk]] <- mdl.filtered$Ct
		n0.k[[kk]] <- mdl.filtered$nt
		s0.k[[kk]] <- mdl.filtered$st
		
		#	Collecting the m_t´s in 'm_models'
		m_models[rownames(mdl.filtered$mt),kk] <- mdl.filtered$mt
		
		#	Collecting information
		y.filt[tt,kk] <- mdl.filtered$ft
		Q.filt[tt,kk] <- mdl.filtered$Qt
		
		}	#	END of "KK"
		
	#
	#	Updating dynamic model probabilities
	#
	if(tt == 1){
		pi0 <- t(as.matrix(rep(1/KK, KK)))			#	Pior model probabilities	
		colnames(pi0) <- c(1:KK)
		rownames(pi0) <- c("Pi_0")
		}else{
		}
	y.filt.v <- y.filt[tt,]
	Q.filt.v <- Q.filt[tt,]
	pit <- dlmbma_postfilter(pi0, gammaa, eps,  yt, y.filt.v, Q.filt.v)
	
	#	Collecting information
	pi.filt[tt,] <- pit
	pi0 <- pit
	rownames(pi0) <- c("pi0")

	# 	Model averaging prediction (updated predictions)
	y.filt.ma[tt] <- pit %*% as.matrix(y.filt[tt,])
	Q.filt.ma[tt] <- pit %*% as.matrix(Q.filt[tt,])
	
	#	Two scenarios
	y.filt.ma.one[tt] <- pit[,c(1:modelthreshold)] %*% as.matrix(y.filt[tt,c(1:modelthreshold)]) / sum(pit[,c(1:modelthreshold)])
	Q.filt.ma.one[tt] <- pit[,c(1:modelthreshold)] %*% as.matrix(Q.filt[tt,c(1:modelthreshold)]) / sum(pit[,c(1:modelthreshold)])
	y.filt.ma.two[tt] <- pit[,c((modelthreshold+1):KK)] %*% as.matrix(y.filt[tt,c((modelthreshold+1):KK)]) / sum(pit[,c((modelthreshold+1):KK)])
	Q.filt.ma.two[tt] <- pit[,c((modelthreshold+1):KK)] %*% as.matrix(Q.filt[tt,c((modelthreshold+1):KK)]) / sum(pit[,c((modelthreshold+1):KK)])
	
	# The m_t´s
	m.filt.ma[tt,] <- t(m_models[,] %*% pit[,] / sum(pit[,]))
	m.filt.ma.one[tt,] <- t(m_models[,c(1:modelthreshold)] %*% pit[,c(1:modelthreshold)] / sum(pit[,c(1:modelthreshold)]))
	m.filt.ma.two[tt,] <- t(m_models[,c((modelthreshold+1):KK)] %*% pit[,c((modelthreshold+1):KK)] / sum(pit[,c((modelthreshold+1):KK)]))
	}	#	END of "TT" 

#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#	Backward Smoothing
#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
tt <- TT-1	
k <- 1
for(tt in (TT-1):1){
	for (k in 1:KK){
		
		}	#	END of "KK"
	} 	#	END of "TT"
#
#	Output
#
list(y.filt = y.filt, Q.filt= Q.filt,
	 pi.filt = pi.filt,
	 y.filt.ma = y.filt.ma, Q.filt.ma = Q.filt.ma,
	 y.filt.ma.one = y.filt.ma.one, Q.filt.ma.one = Q.filt.ma.one,
	 y.filt.ma.two = y.filt.ma.two, Q.filt.ma.two = Q.filt.ma.two,
	 m.filt.ma = m.filt.ma, 
	 m.filt.ma.one = m.filt.ma.one,
	 m.filt.ma.two = m.filt.ma.two
	)
}
#
#	END of "dlmbma.R"