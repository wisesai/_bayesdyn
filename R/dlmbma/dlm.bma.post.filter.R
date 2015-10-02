dlm.bma.post.filter <- function(pi0, gammaa, eps,  yt, y.filt.v, Q.filt.v){
#
#	Computing dynamic posterior model probabilities
#
#	Inputs:
#  	pi0			-	K-vector of input model probabilities
#	gammaa		-	flattening parameter
#	eps			-	minimum threshold for model probabilities  
#	yt			-	observed value of y_t
#	y.filt.v	-	K-vector of predicted values of y_t | y_{t-1} from the Kalman filter
#	Q.filt.v	-	K-vector of predicted variances of y_t | y_{t-1} from the Kalman filter
#
# 	Output:
# 	pit			-	K-vector of updated model probabilities
#
#	Author: 	Juan~Carlos Martínez-Ovando
#	Email:		juan.martinez@banxico.org.mx
#				JC.Martinez.Ovando@gmail.com
#
#	Reference:	Martínez-Ovando, J.~C. (2014) "The Influence of Crime on the Economic Growth 
#					Dynamics in Mexico in the Short-Term," Banco de México, Mimeo.
#

# Form predicted pi values
pit <- pi0^gammaa / sum(pi0^gammaa)

# Update pi values
logpyt <- -0.5*log(Q.filt.v) - 0.5*(yt-y.filt.v)^2/Q.filt.v
logpyt <- logpyt - max(logpyt)
pyt <- exp (logpyt)
pit <- Q.filt.v * pyt
pit <- pit/sum(pit)
pit <- pit + eps
pit <- pit/sum(pit)
pit <- t(as.matrix(pit))
rownames(pit) <- c("pit")

# Output
return(pit)
}
#
#	END of "dlm.bma.post.filter.R"