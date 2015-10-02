ar <- function(x,p){
#
#	This function fits a reference AR(p) model to a time-series data
#
#	Reference:	Prado & West (2010)
#
#	Inputs:
#   	x		-	Time series data with dimension (Tx1)
#  		p		-	Model order
#
#	Outputs: 
#  		m  		-	...
#  		r		-	...
#  		nu		-	...
#  		s		-	...
#

#	Organizing the data
arx <- x-mean(x)
mi <- min(arx)
ma <- max(arx) 
T <- length(arx)

y <- as.matrix(arx[(p+1):T])
X <- matrix(NA,nrow =(T-p), ncol=p)
for(t in T:(p+1)){X[t-p,] <- arx[(t-1):(t-p)]}

#	Estimation
xtx <- t(X)%*%X
m <- solve(xtx) %*% t(X)%*%y
r <- y - X%*%m
nu <- T - 2*p
s <- solve(nu)*(t(r)%*%r)

# 	Output
output <- list(m,r,nu,s)
return(output)
#
#	--	END	--
}