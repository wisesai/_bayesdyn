gibbs_A <- function(M,t,y,Z,Btdraw,sigt,Sblockdraw,A_0_prmean,A_0_prvar,numa,S_prmean,S_prvar){
#
#	gibbs_A.R
#
#	Author: 	Juan~Carlos Martínez-Ovando
#	Email:		juan.martinez@banxico.org.mx
#				jc.martinez.ovando@gmail.com
#

#	Substract from the data y(t), the mean Z x B(t)
yhat <- matrix(0, M, t)
i <- 1
for(i in 1:t){
    yhat[ ,i] <- as.matrix(y[ ,i]) - Z[((i-1)*M+1):(i*M),] %*% as.matrix(Btdraw[ ,i])
  }

#	This part is more tricky, check Primiceri, Zc is a [M x M(M-1)/2] matrix
Zc <- - t(yhat[,])
sigma2temp <- matrix(0,t,dim(sigt)[2])
i <- 1
for(i in 1:t){
    sigmatemp <- t(diag(sigt[((i-1)*M+1):(i*M), ]^2))
    sigma2temp[i, ] <- sigmatemp
	}  
    
ind <- 1
ii <- 2
for(ii in 2:M){
	#	Draw each block of A(t)
	y2 <- t(as.matrix(yhat[ii,]))
	Z2 <- as.matrix(Zc[ ,1:(ii-1)])
	Ht2 <- as.matrix(sigma2temp[ ,ii])
	Qt2 <- as.matrix(Sblockdraw[[ii-1]])
	m2 <- sizeS[ii-1]
	p2 <- 1
	t2 <- t
	B02 <- as.matrix(A_0_prmean[((ii-1)+(ii-3)*(ii-2)/2):ind, ])
	V02 <- as.matrix(A_0_prvar[((ii-1)+(ii-3)*(ii-2)/2):ind,((ii-1)+(ii-3)*(ii-2)/2):ind])
									
	gibbs_ssb_out <- gibbs_ssb(y2, Z2, Ht2, Qt2, m2, p2, t2, B02, V02)
    Atblockdraw <- gibbs_ssb_out[[1]]
	log_lik2a <- gibbs_ssb_out[[2]] 
	
	if(ii == 2){
		Atdraw <- Atblockdraw							#	Atdraw is the final matrix of draws of A(t)
	}else{
		Atdraw <- rbind(Atdraw, Atblockdraw)
		}
    ind <- ind + ii
  }
    
#	Sampling S -- the covariance of A(t) (from iWishart).
#	Take the SSE in the state equation of A(t).
Attemp <- t(Atdraw[ ,2:t]) - t(Atdraw[ ,1:(t-1)])
sse_2 <- matrix(0, numa, numa)
i <- 1
for(i in 1:(t-1)){
    sse_2 <- sse_2 + as.matrix(Attemp[i,]) %*% t(as.matrix(Attemp[i, ]))
  }

#	...and subsequently draw S, the covariance matrix of A(t) 
ijc <- 1
jj <- 2
for(jj in 2:M){
    Sinv <- solve(sse_2[((jj-1)+(jj-3)*(jj-2)/2):ijc, ((jj-1)+(jj-3)*(jj-2)/2):ijc] + S_prmean[[jj-1]]);
    Sinvblockdraw <- wishartrnd(Sinv,t+S_prvar[jj-1])
    Sblockdraw[[jj-1]] <- solve(Sinvblockdraw)			#	Sblockdraw this is a draw from S
    ijc <- ijc + jj
  }

#	Output
output <- list(Atdraw,Sinvblockdraw,yhat)
return(output)
}
#
#	-- End of "gibbs_A.R" --