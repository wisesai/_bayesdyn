tvar.pred <- function(x,p,K,Msim,m,C,n,s){
#
#	This function generates predictions from a TVAR(p) fitted model
#	"k" periods ahead, using direct simulations methods
#
#	Reference:	West & Harrison (1997), Prado & West (1997)
#
#	Inputs:
#   	x		-	Time series data with dimension (Tx1)
#  		p		-	Model order
#  		K		-	Periods ahead for prediction
#  		Msim	-	number of simulations from the predictive distribution
#  		m  		-	(pxT) matrix with posterior mean vector states
#  		C		-	(pxpxT) post variance matrix for states
#  		n		-	T posterior df´s for observation variances
#  		s		-	T posterior observation variances estimates
#
#	Outputs: 
#  		x_Kpred	- 	(KxM) matrix with predictive simulations from 
#					the jint predictive distribution of "K" periods ahead
#

#	Organizing the data
arx <- x-mean(x)
T <- nrow(arx)
x_Kpred <- matrix(NA,(T+K),Msim)

#	Inicio de los ajustes (variables de estados)
mT <- m[T,]
CT <- C[T,,]
sT <- s[T] 
nT <- n[T]

#	Prediction
msim <- 1
for(msim in 1:Msim){
	#	Simulating state-variables
	arxp <- rbind(arx,matrix(NA,nrow=K,ncol=1))
	mTfut <- mvrnorm(n=1, mT, CT)
	#mTfut <- mT
	sTfut <- rgamma(n=1, shape=sT, rate = nT)
	#	Random trayectories
	for(k in 1:K){
		F <- as.matrix(arxp[(T+k-1):(T+k-1-p+1),])
		arxp[T+k,1] <- mvrnorm(n=1, mTfut%*%F, sTfut)
		}
	x_Kpred[,msim] <- arxp 
	}
x_Kpred <- ts(x_Kpred,frequency = 12, start = c(2003, 2))

# 	Output
output <- x_Kpred
return(output)
#
#	--	END	--
}