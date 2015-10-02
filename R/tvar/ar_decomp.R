ar_decomp <- function(x,p,m){
#
#	This function decomposes the parameters of an AR(p) model in webcomponents
#
#	Reference:	Prado & West (2010)
#
#	Inputs:
#   	x		-	Time series data with dimension (Tx1)
#  		p		-	Model order
#		m 		-	AR(p) coeficient parameter, with dimension (1xp)
#
#	Outputs: 
#   	waves	- the estimated wavelengths of ARMA(2,1) components with complex roots
#   	mods	- the estimated moduli of all components, the first 
#   			  ones corresponding to complex roots in order in waves	
#   	decomp	- the matrix of estimated components, in same order from top town, 
#                 starting at t=p and running to t=length(x) 
#

#	Organizing the data
arx <- x-mean(x)
mi <- min(arx)
ma <- max(arx) 
T <- length(arx)

y <- as.matrix(arx[(p+1):T])
X <- matrix(NA,nrow =(T-p), ncol=p)
for(t in T:(p+1)){X[t-p,] <- arx[(t-1):(t-p)]}

# 	Decomposition ... (this)
G <- matrix(0,p,p)
G[1,] <- m
G[2:p,1:(p-1)] <- diag(p-1)
eigG <- eigen(G)

#	#	Pendiente
#	[e,l]=eig(G); l=diag(l); arg=angle(l); mods=abs(l); H=diag(e(1,:))*inv(e);
#

# 	Output
output <- list(b,r,nu,s)
return(output)
#
#	--	END	--
}