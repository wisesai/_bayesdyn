gibbs_logSigma <- function(M,t,Atdraw,yhat,statedraw,u2_s,m_s,Zs,Wdraw,sigma_prmean,sigma_prvar,q_s,W_prmean,W_prvar){
#
#	gibbs_logSigma.R
#
#	Author: 	Juan~Carlos Martínez-Ovando
#	Email:		juan.martinez@banxico.org.mx
#				jc.martinez.ovando@gmail.com
#

#	First create capAt, the lower-triangular matrix A(t) with ones on the main diagonal. 
#	This is because the vector Atdraw has a draw of only the non-zero and non-one elements of A(t) only.
capAt <- matrix(0, M*t, M)
i <- 1
for(i in 1:t){
    capatemp <- diag(M)
    aatemp <- as.matrix(Atdraw[ ,i])
    ic <- 1
	j <- 2
    for(j in 2:M){
        capatemp[j, 1:(j-1)] <- t(aatemp[ic:(ic+j-2),1]);
        ic <- ic + j - 1
	  }#	-- End of "for 2"
	capAt[((i-1)*M+1):(i*M), ] <- capatemp
  }#	-- End of "for 1"


#	yhat is the vector y(t) - Z x B(t) defined previously. Multiply yhat with capAt, 
#	i.e the lower triangular matrix A(t). Then take squares of the resulting quantity (saved in matrix y2)
i <- 1
for(i in 1:t){
    ytemps <- capAt[((i-1)*M+1):(i*M), ] %*% yhat[ ,i]
    if(i == 1){
		y2 <- ytemps^2
	  }else{
		y2 <- cbind(y2,(ytemps^2))
	  }
	}
	
#	Transform to the log-scale but also add the 'offset constant' to prevent the case where y2 
#	is zero (in this case, its log would be -Infinity) 
yss <- matrix(0, t, M)
i <- 1
for(i in 1:M){
    yss[ ,i] <- log(t(y2[i,]) + 0.001)
  }
	
#	In order to draw the log-volatilies, substract the mean and variance of the 7-component 
#	mixture of Normal approximation to the measurement error covariance
vart <- matrix(0, t*M, M)
yss1 <- matrix(0, t, M)
i <- 1
for(i in 1:t){
	j <- 1
    for(j in 1:M){
        imix <- statedraw[i, j]
        vart[(i-1)*M+j,j] <- u2_s[imix]
        yss1[i,j] <- yss[i, j] - m_s[imix] + 1.2704
	  } # -- End of "for 2"
  } # -- End of "for 1"

#	Sigtdraw is a draw of the diagonal log-volatilies, which will give us Sigma(t)
y3 <- t(yss1)
Z3 <- Zs
Ht3 <- vart
Qt3 <- Wdraw
m3 <- M
p3 <- M
t3 <- t
B03 <- sigma_prmean
V03 <- sigma_prvar
gibbs_ss_out3 <- gibbs_ss(y3, Z3, Ht3, Qt3, m3, p3, t3, B03, V03);
Sigtdraw <- gibbs_ss_out3[[1]]
log_lik3 <- gibbs_ss_out3[[2]]	

#	Next draw statedraw (chi square approximation mixture component) conditional on Sigtdraw
#	This is used to update at the next step the log-volatilities Sigtdraw
jj <- 1
for(jj in 1:M){
	i <- 1
    for(i in 1:t){
		j <- 1
        for(j in 1:length(m_s)){
            temp1 <- (1/sqrt(2*pi*u2_s[j]))*exp(-.5*(((yss[i,jj] - Sigtdraw[jj,i] - m_s[j] + 1.2704)^2)/u2_s[j]));
            prw[j, 1] <- q_s[j,1]*temp1;
		  } # -- End of "for 3"
        prw <- prw / sum(prw)
        cprw <- as.matrix(cumsum(prw))
        trand <- runif(1,0,1)								#	rand(1,1); % in Matlab 
        if(trand < cprw[1,1]){
			imix <- 1
		  }else if(trand < cprw[2,1]){
			imix <- 2
		  }else if(trand < cprw[3,1]){
			imix<- 3
		  }else if(trand < cprw[4,1]){
			imix <- 4
		  }else if(trand < cprw[5,1]){
			imix <- 5
		  }else if(trand < cprw[6,1]){
			imix <- 6
		  }else{
			imix <- 7
		  }
        statedraw[i,jj] <- imix								#	this is a draw of the mixture component index
	  }# -- End of "for 2"
  } # -- End of "for 1"

#	Draws in Sigtdraw are in logarithmic scale (log-volatilies). 
#	Create original standard deviations of the VAR covariance matrix.
sigtemp <- diag(M)
sigt <- matrix(0, M*t, M)
i <- 1
for(i in 1:t){
	j <- 1
    for(j in 1:M){
        sigtemp[j,j] <- exp(0.5*Sigtdraw[j,i])
	  } # -- End of "for 2"
	sigt[((i-1)*M+1):(i*M), ] <- sigtemp
  } # -- End of "for 1"

#	Sampling W -- the covariance of SIGMA(t) (from iWishart).
#	Get first differences of Sigtdraw to compute the SSE.
Sigttemp <- t(Sigtdraw[,2:t]) - t(Sigtdraw[ ,1:(t-1)])

sse_2 <- matrix(0, M, M)
i <- 1
for(i in 1:(t-1)){
    sse_2 <- sse_2 + Sigttemp[i, ] %*% t(Sigttemp[i, ])
  }
Winv <- solve(sse_2 + W_prmean)
Winvdraw <- wishartrnd(Winv, t + W_prvar)
Wdraw <- solve(Winvdraw)									#	this is a draw from W		
	
#	Create the VAR covariance matrix H(t). It holds that:
#           A(t) x H(t) x A(t)' = SIGMA(t) x SIGMA(t) '
Ht <- matrix(0, M*t, M)
Htsd <- matrix(0, M*t, M)
i <- 1
for(i in 1:t){
	inva <- solve(capAt[((i-1)*M+1):(i*M), ])
    stem <- sigt[((i-1)*M+1):(i*M), ]
    Hsd <- inva %*% stem
    Hdraw <- Hsd %*% t(Hsd)
    Ht[((i-1)*M+1):(i*M), ] <- Hdraw						#	this is a draw from H(t)
    Htsd[((i-1)*M+1):(i*M), ] <- Hsd						#	Cholesky of H(t)
  }

#	Output
output <- list(Wdraw,Ht,Htsd)
return(output)
}
#
#	-- End of "gibbs_logSigma.R" --