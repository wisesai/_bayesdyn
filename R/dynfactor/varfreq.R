varfreq <- function(ydata,lags,xdata){
#
#	varfreq.R
#	This function computes frequentist estimates for VAR models
#
#	Original code was written by Chris Sims and adapted by J.C. Martinez-Ovando (2015)
#
#	Author: 	Juan~Carlos Martínez-Ovando
#	Email:		juan.martinez@banxico.org.mx
#				jc.martinez.ovando@gmail.com
#

if(lags != 1){
  stop('estvar: Currently available only for lags==1.')
  }

t <- dim(ydata)[1]
ny <- dim(ydata)[2]
nox <- is.null(xdata)
if(nox == FALSE){
    t1 <- dim(xdata)[1]
    nx <- dim(xdata)[2]
  }else{
   t1 <- t
   nx <- 0
  }

if(t1 != t){
  print('estvar: Mismatch of x and y data lengths')
  }

if(lags == 1){
	X <- as.matrix(ydata[(lags+1-1):(t-1), ])
#	X <- matrix(0, c(t, ny) ) 						#	array(0, c(t, ny, lags) ) 					# 	in Matlab for multiple 'lags'
#	i <- 1
#	for(i in 1:lags){
#		X[(lags+1):t, ] <- ydata[(lags+1-i):(t-i), ]
#		}
	if(nox == FALSE){
		X <- cbind(X, as.matrix(xdata[(lags+1):t, ]))
	  }else{
		X <- X
	  }
	y <- as.matrix(ydata[(lags+1):t, ])
#	X <- X[(lags+1):t, ]

	k <- ny*lags+nx									#	num. of coefficients
	svd_out <- svd(X)								#	singular value decomposition, X is t*k
	if(length(svd_out$d) > 1){
		d <- diag(svd_out$d)
	  }else if(length(svd_out$d) == 1){
		d <- matrix(svd_out$d,1,1)
	  }
	vr <- svd_out$v 
	vl <- svd_out$u
	
	if(length(svd_out$d) > 1){
		d <- diag(1 / diag(d))
	  }else if(length(svd_out$d) == 1){
		d <- 1 / d
	  }
	B <- t(vl) %*% y
	B <- (vr * (matrix(1,k,dim(d)[2])%*%d)) %*% B	#	this is OLS equ by equ, B is k*ny
	u <- y - X %*% B
	omega <- (t(u)%*%u) * (t-lags)^-1				#	covarinace matrix
	xx <- vr * (matrix(1,k,dim(d)[2])%*%d)
	xx <- xx %*% t(xx)
	By <- B[1:(ny*lags), ]
#	By.temp <- array(0, c(ny,lags,ny)) 
#	i <- 1
#	for(i in i:ny){									#	By <- reshaper(By,'col')						#	variables, lags, equation:
#		By.temp[, , i] <- By[, i]
#		}
	By <- t(By)										#	equations, variables, lags
  }
  
if(nox == TRUE){
	Bx <- NULL
  }else{
	Bx <- t(B[(ny*lags+(1:nx)), ])					#	equations*nx
  }

output <- list(By,Bx,u,omega,xx)
return(output)
}
#
# -- End of "varfreq.R"