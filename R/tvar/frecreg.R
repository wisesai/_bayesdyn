frecreg <- function(y,X){
#
#	This function fits a Bayesian regression model
#	with a conjugate informative prior
#

T <- nrow(X)
p <- ncol(X)
 
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