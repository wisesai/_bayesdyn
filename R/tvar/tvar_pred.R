tvar <- function(x,p,del,m0,C0,s0,n0){
#
#	This function fits a TVAR(p) model written in DLM form,
#	with discount factors "del(1)" for state and "del(2)" for obs var
#
#	Reference:	West & Harrison (1997), Prado & West (1997)
#
#	Inputs:
#   	x		-	Time series data with dimension (Tx1)
#  		p		-	Model order
#  		m0 		-	(px1) vector prior mean for state variables
#  		C0 		-	(pxp) prior variance matrix for state variables
#  		n0		- 	Prior df for observation variance
#  		s0 		-	Prior estimate of observation variance
#
#	Outputs: 
#  		m  		-	(pxT) matrix with posterior mean vector states
#  		C		-	(pxpxT) post variance matrix for states
#  		n		-	T posterior df´s for observation variances
#  		s		-	T posterior observation variances estimates
#  		e  		-	(Tx1) estimated innovations (zero up to t=p) 
#

#	Organizing the data
d <- del[1]
b <- del[2];
arx <- x-mean(x)
T <- nrow(arx)

#	Initial states for model
m <- matrix(NA,nrow=p,ncol=p)
rownames(m) <- c(1:p)
colnames(m) <- c(1:p)
for(t in 1:p){m[t,] <- m0}
C <- array(NA, dim=c(T,p,p))
for(t in 1:p){C[t,,] <- C0}
s <- s0 * matrix(1,nrow=p,ncol=1)
n <- n0 * matrix(1,nrow=p,ncol=1)

#	Inicio de los ajustes (variables de estados)
mt <- m0
Ct <- C0
st <- s0 
nt <- n0

# 	Forward filtering
t <- p+1
for(t in (p+1):T){
	F <- as.matrix(arx[(t-1):(t-p),])
    A <- Ct%*%F/d
	q <- t(F)%*%A + st
	A <- as.numeric(solve(q))*A
	e <- arx[t]-t(F)%*%mt	
	
	#	Updating state variables
	mt <- mt + A%*%e
	m <- rbind(m,t(mt))
    r <- b*nt+e*e/q
	nt <- b*nt+1
	r <- r/nt
	st <- st*r
	n[t] <- nt
	s[t]=st
	Ct <- as.numeric(r)*(Ct/d-(A%*%t(A))*as.numeric(q))
	Ct <- (Ct+ t(Ct))/2
    C[t,,] <- Ct
	}

#	Backward smoothing
e <- matrix(0,T,1)
t <- (T-1)
for(t in (T-1):(p+1)){
    m[t,] <- (1-d)*m[t,] + d*m[t+1,]
    if(t > p){
        e[t] = arx[t]-m[t,]%*%arx[(t-1):(t-p),]
    }
    n[t] <- (1-b)*n[t]+ b*n[t+1]
	st <- s[t]
	s[t] <- 1/((1-b)/st+b/s[t+1])
	C[t,,] <- s[t]*((1-d)*C[t,,]/st + d*d*C[t+1,,]/s[t+1])
	}

# 	Ad-hoc treatment of first "p" values 
for(t in 1:p){m[t,] <- m[p+1,]; C[t,,] <- C[p+1,,]; n[t] <- n[p+1]; s[t] <- s[p+1]}

# 	Output
output <- list(m,C,s,n,e)
return(output)
#
#	--	END	--
}