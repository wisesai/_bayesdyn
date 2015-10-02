bvar <- function(Yraw, constant, p_lags, forecasting, forecast_method, N_pred_h, impulses, impulses_h, prior, N_mcmc, N_burn, N_print){
#
#	bvar.R
#	VAR using the Gibbs sampler, based on independent Normal Wishar prior
#	--------------------------------------------------------------------------
#	Bayesian estimation, prediction and impulse response analysis in VAR
#	models using the Gibbs sampler. Dependent on your choice of forecasting,
#	the VAR model is:
#
#	Iterated forecasts:
#   	  Y(t) = A0 + Y(t-1) x A1 + ... + Y(t-p_lags) x Ap + e(t)
#
#	so that in this case there are p_lags lags of Y (from 1 to p_lags).
#
#	Direct N_pred_h-step ahead foreacsts:
#   	  Y(t+N_pred_h) = A0 + Y(t) x A1 + ... + Y(t-p_lags+1) x Ap + e(t+h)
#
#	so that in this case there are also p_lags lags of Y (from 0 to p_lags-1).
#
#	In any of the two cases, the model is written as:
#
#                   Y(t) = X(t) x A + e(t)
#
#	where e(t) ~ N(0,SIGMA), and A summarizes all parameters. Note that we
#	also use the vector a which is defined as a=vec(A).
#
#	Original code was written by D. Korobilis and adpated by J.C. Martinez-Ovando (2015)
#
#	Author: 	Juan~Carlos Martínez-Ovando
#	Email:		juan.martinez@banxico.org.mx
#				jc.martinez.ovando@gmail.com
#

#	--------------------------  Setting up data  -----------------------------------
#	Get initial dimensions of dependent variable
Traw <- dim(Yraw)[1]
M <- dim(Yraw)[2]

#	The model specification is different when implementing direct forecasts,
#	compared to the specification when computing iterated forecasts.
if(forecasting == 1){
    if(N_pred_h <= 0){
        print('You have set forecasting, but the forecast horizon choice is wrong')
      }

    #	Now create VAR specification according to forecast method
    if(forecast_method == 0){				#	Direct forecasts
        Y1 <- Yraw[N_pred_h+1:dim(Yraw)[1], ]
        Y2 <- Yraw[2:(dim(Yraw)[1]-N_pred_h), ]
        Traw <- Traw - N_pred_h - 1
	  }else if(forecast_method == 1){		#	Iterated forecasts
        Y1 <- Yraw
        Y2 <- Yraw
      }else{
        print('Wrong choice of forecast_method')
      }
	}else{
		Y1 <- Yraw
		Y2 <- Yraw
	}

#	Generate lagged Y matrix. This will be part of the X matrix
Ylag <- mlag2(Y2,p_lags)							#	Y is [T x M]. ylag is [T x (Mp)]

#	Now define matrix X which has all the R.H.S. variables (constant, lags of
#	the dependent variable and exogenous regressors/dummies)
if(constant == 1){
    X1 <- cbind( matrix(1, Traw-p_lags,1), Ylag[(p_lags+1):Traw, ])
  }else{
    X1 <-  Ylag[(p_lags+1):Traw, ]
  }

#	Get size of final matrix X
Traw3 <- size(X1)[1] 
K <- size(X1)[2]

#	Create the block diagonal matrix Z
Z1 <- kronecker(diag(M), X1)

#	Form Y matrix accordingly
#	Delete first "LAGS" rows to match the dimensions of X matrix
Y1 <- Y1[(p_lags+1):Traw, ]						#	This is the final Y matrix used for the VAR

#	Traw was the dimesnion of the initial data. T is the number of actual 
#	time series observations of Y and X (we lose the p_lags-lags)
T <- Traw - p_lags

#	-----------------------  Setting up forecasting  --------------------------------
#	Now keep also the last "N_pred_h" or 1 observations to evaluate (pseudo-)forecasts
if(forecasting == 1){
    Y_pred <- matrix(0, N_mcmc*N_pred_thin, M)	#	Matrix to save prediction draws
    PL <- matrix(0, N_mcmc, 1)				#	Matrix to save Predictive Likelihood
    
    if(forecast_method == 0){				#	Direct forecasts, we only need to keep the 
        Y <- Y1[1:(dim(Y1)[1]-1),	]		#	last observation
        X <- X1[1:(dim(X1)[1]-1),	]
        Z <- kronecker(diag(M), X)
        T <- T - 1
	  }else{								#	Iterated forecasts, we keep the last N_pred_h observations
        Y <- Y1[1:(dim(Y1)[1]-N_pred_h), ]
        X <- X1[1:(dim(X1)[1]-N_pred_h), ]
        Z <- kronecker(diag(M), X)
        T <- T - N_pred_h;
      }
	}else{
    Y <- Y1
    X <- X1
    Z <- Z1
	}

#	------------------  Setting up impulse-responses  ---------------------------
#	Creating matrices to store forecasts
if(impulses == 1){
    #	Impulse response horizon
    imp <- array(0, c(N_mcmc, M, impulses_h))
    bigj <- matrix(0, M, M*p_lags)
    bigj[1:M, 1:M] <- diag(M)
	}


#	----------------------------  Preliinaries  --------------------------------
#	First get ML estimators
Y <- as.matrix(Y)
X <- as.matrix(X)
A_OLS <- solve(t(X)%*%X) %*% t(X) %*% Y			#	This is the matrix of regression coefficients
a_OLS <- A_OLS									#	A_OLS(:);         % This is the vector of parameters, i.e. it holds
												#	that a_OLS = vec(A_OLS)
SSE <- t(Y - X %*% A_OLS) %*% (Y - X %*% A_OLS)	#	Sum of squared errors
SIGMA_OLS <- SSE / (T-K+1)

#	Initialize Bayesian posterior parameters using OLS values
alpha <- a_OLS									#	This is the single draw from the posterior of alpha
ALPHA <- A_OLS									#	This is the single draw from the posterior of ALPHA
SSE_Gibbs <- SSE								#	This is the single draw from the posterior of SSE
SIGMA <- SIGMA_OLS								#	This is the single draw from the posterior of SIGMA

#	Storage space for posterior draws
alpha_draws <- matrix(0, N_mcmc, K*M)
ALPHA_draws <- array(0, c(N_mcmc, K, M))
SIGMA_draws <- array(0, c(N_mcmc,M,M))

#	Storage space for posterior draws
alpha_draws <- matrix(0, N_mcmc, K*M)
ALPHA_draws <- array(0, c(N_mcmc,K,M))
SIGMA_draws <- array(0, c(N_mcmc,M,M))

#	-----------------  Prior specification  --------------------------------------
n <- K*M										#	Total number of parameters (size of vector alpha)
#	Define hyperparameters
if(prior == 1){									#	Normal-Wishart
    a_prior <- matrix(0, n, 1)					#	<---- prior mean of alpha (parameter vector)
    V_prior <- 10*diag(n)						#	<---- prior variance of alpha
    
    #	Hyperparameters on inv(SIGMA) ~ W(v_prior,inv(S_prior))
    v_prior <- M								#	<---- prior Degrees of Freedom (DoF) of SIGMA
    S_prior <- diag(M)							#	<---- prior scale of SIGMA
    inv_S_prior <- solve(S_prior)
  }else if(prior == 2 ){						#	Minnesota-Whishart
    #	Prior mean on VAR regression coefficients
    A_prior <- rbind(matrix(0, 1, M), 0.9*diag(M), matrix(0, ((p_lags-1)*M), M))	#	<---- prior mean of ALPHA (parameter matrix) 
    a_prior <- matrixcol(A_prior)											#	<---- prior mean of alpha (parameter vector)
    
    #	Minnesota Variance on VAR regression coefficients
    #	First define the hyperparameters 'a_bar_i'
    a_bar_1 <- 0.5
    a_bar_2 <- 0.5
    a_bar_3 <- 10^2
    #	Now get residual variances of univariate p_lags-lag autoregressions. Here
    #	we just run the AR(p_lags) model on each equation, ignoring the constant
    #	and exogenous variables (if they have been specified for the original
    #	VAR model)
    sigma_sq <- matrix(0, M, 1)					#	vector to store residual variances
	i <- 1
    for(i in 1:M){
        #	Create lags of dependent variable in i-th equation
        Ylag_i <- mlag2(Yraw[,i],p_lags)
        Ylag_i <- Ylag_i[(p_lags+1):Traw, ]
        #	Dependent variable in i-th equation
        Y_i <- Yraw[(p_lags+1):Traw, i]
        #	OLS estimates of i-th equation
        alpha_i <- solve(t(Ylag_i)%*%Ylag_i)%*%(t(Ylag_i)%*%Y_i)
        sigma_sq[i, 1] <- (1 /(Traw-p_lags+1))*t(Y_i - Ylag_i%*%alpha_i)%*%(Y_i - Ylag_i%*%alpha_i);
		} # -- End of "for"

    #	Now define prior hyperparameters.
    #	Create an array of dimensions K x M, which will contain the K diagonal
    #	elements of the covariance matrix, in each of the M equations.
    V_i <- matrix(0, K, M)
	#	index in each equation which are the own lags
    ind <- matrix(0, M, p_lags)
	i <- 1
    for(i in 1:M){
        ind[i, ] <- indexownlags(i,constant,M,p_lags,K)
		}
    i <- 1;
	for(i in 1:M){														#	for 1 - each i-th equation
        j <- 1;
		for(j in 1:K){													#	for 2 - each j-th RHS variable
            if(constant == 1){#	if there is constant in the model use this code:
                A <- ind[i, ] 
				if(j == 1){
                    V_i[j, i] <- a_bar_3*sigma_sq[i,1]					#	variance on constant                
                }else if( finding(ind[i, ],j) > 0 ){
                    V_i[j, i] <- a_bar_1 / (p_lags^2)					#	variance on own lags           
                }else{
					kj <- 1
                    for(kj in 1:M){
                        if( finding(ind[kj, ],j) > 0 ){
                            ll <- kj                   
							}
						}												#	variance on other lags   
                    V_i[j, i] <- (a_bar_2*sigma_sq[i, 1]) / ((p_lags^2)*sigma_sq[ll,1])      
                }
            }else{														#	if there is no constant in the model use this:
                if( finding(ind[i, ],j) >0 ){
                    V_i[j, i] <- a_bar_1 / (p_lags^2)					#	variance on own lags
                }else{
					kj <- 1
                    for(kj in 1:M){
                        if( finding(ind[kj, ],j) > 0 ){
                            ll <- kj
							}                        
						}# variance on other lags  
                    V_i[j, i] <- (a_bar_2*sigma_sq[i,1]) / ((p_lags^2)*sigma_sq[ll,1])            
					}
				}
			}# -- End of "for 2"
		}# -- End of "for 1"

	#	Now V is a diagonal matrix with diagonal elements the V_i
    V_prior <- matrix(0,length(matrixcol(V_i))[1],length(matrixcol(V_i))[1])
	diag(V_prior) <- matrixcol(V_i)										#	this is the prior variance of the vector alpha

    #	Hyperparameters on inv(SIGMA) ~ W(v_prior,inv(S_prior))
    v_prior <- M
    S_prior <- diag(M)
    inv_S_prior <- solve(S_prior)   
  }# --	End of "else if"
  
#	-----------------  MCMC sampling  --------------------------------------
tic <- proc.time()														#	This is just a timer
print('Number of iterations')
mod <- function(x,m){t1<-floor(x/m); return(x-t1*m)}

irep <- 1
for(irep in 1:(N_mcmc + N_burn)){										#	GIBBS iterations starts here
	#	Print iterations
    if( (mod(irep,N_print) == 0) & (irep < N_burn) ){
	    print(paste("Burn-in: ", irep, " of ", N_burn, sep= ""))
	  }else if( (mod(irep,N_print) == 0) & (irep > N_burn) ){
	    print(paste("MCMC: ", irep-N_burn, " of ", N_mcmc, sep= ""))
	  }

    VARIANCE <- kronecker(solve(SIGMA),diag(T))
    V_post <- solve(solve(V_prior) + t(Z) %*% VARIANCE %*% Z)
    a_post <- V_post %*% (solve(V_prior) %*% a_prior + t(Z) %*% VARIANCE %*% matrixcol(Y))
    alpha <- a_post + t(chol(V_post)) %*% as.matrix(rnorm(n, mean=0, sd=1))		#	Draw of alpha
	
    ALPHA <- reshapes(alpha,'col',K,M)											#	reshape(alpha,K,M) # in Matlab  -- Draw of ALPHA
    
    #	Posterior of SIGMA|ALPHA,Data ~ iW(inv(S_post),v_post)
    v_post <- T + v_prior
    S_post <- S_prior + t(Y - X%*%ALPHA) %*% (Y - X%*%ALPHA)
	SIGMA <- solve(wishartrnd(solve(S_post), v_post))									#	Draw SIGMA

    #	Store results  
    if(irep > N_burn){              
        #	=========FORECASTING:
        if(forecasting == 1){
            if(forecast_method == 0){											#	Direct forecasts
                Y_temp <- zeros(N_pred_thin,M)
                #	compute 'N_pred_thin' predictions for each draw of ALPHA and SIGMA
                ii <- 1
				for(ii in 1:N_pred_thin){
                    X_fore <- cbind(as.matrix(1), t(as.matrix(Y[T, ])), t(as.matrix(X[T, 2:(M*(p_lags-1)+1)])))
                    #	Forecast of T+1 conditional on data at time T
                    Y_temp[ii, ] <- X_fore %*% ALPHA + t(as.matrix(rnorm(M, mean=0, sd=1))) %*% chol(SIGMA)
					}
                #	Matrix of predictions
                Y_pred[ (((irep-N_burn)-1)*N_pred_thin+1):((irep-N_burn)*N_pred_thin), ] <- Y_temp
                #	Predictive likelihood
                PL[(irep-N_burn), ] <- dmvnorm( Y1[T+1, ], X[T, ]%*%ALPHA, SIGMA)
                if(PL[(irep-N_burn), ] == 0){
                    PL[(irep-N_burn), ] = 1;
					}
            }else if(forecast_method == 1){										#	Iterated forecasts
                #	 The usual way is to write the VAR(p_lags) model in companion
                #	form, i.e. as VAR(1) model in order to estimate the
                #	N_pred_h-step ahead forecasts directly (this is similar to the 
                #	code we use below to obtain impulse responses). Here we 
                #	just iterate over N_pred_h periods, obtaining the forecasts at 
                #	T+1, T+2, ..., T+N_pred_h iteratively.
                Y_temp2 <- matrix(NaN,N_pred_thin,M)
				ii <- 1
                for(ii in 1:N_pred_thin){ # -- "for 1"
                    #	Forecast of T+1 conditional on data at time T
                    X_fore <- cbind(as.matrix(1), t(as.matrix(Y[T, ])), t(as.matrix(X[T,2:(M*(p_lags-1)+1)])) )
                    Y_hat <- X_fore%*%ALPHA + t(as.matrix(rnorm(M, mean=0, sd=1))) %*% chol(SIGMA)
                    Y_temp <- Y_hat
                    X_new_temp <- X_fore
					i <- 1
                    for(i in 1:(N_pred_h-1)){ # -- "for 2"								#	Predict T+2, T+3 until T+N_pred_h                   
                        if(i <= p_lags){
                            #	Create matrix of dependent variables for                       
                            #	predictions. Dependent on the horizon, use the previous                       
                            #	forecasts to create the new right-hand side variables
                            #	which is used to evaluate the next forecast.                       
							if(2 < ((M*(p_lags-i))+1)){
                                X_new_temp <- cbind(as.matrix(1),  as.matrix(Y_hat), t(as.matrix(X_fore[,2:((M*(p_lags-i))+1)])) )
							  }else{
							    X_new_temp <- cbind(as.matrix(1),  as.matrix(Y_hat) )
							  }
                            #	This gives the forecast T+i for i=1,..,p_lags                       
                            Y_temp <- X_new_temp%*%ALPHA + t(as.matrix(rnorm(M, mean=0, sd=1))) %*% chol(SIGMA)                     
                            Y_hat <- cbind(Y_hat, Y_temp)
						  }else{
                            X_new_temp <- cbind( as.matrix(1), t(as.matrix(Y_hat[,1:(M*p_lags)])) )
                            Y_temp <- X_new_temp%*%ALPHA + t(as.matrix(rnorm(M, mean=0, sd=1))) %*% chol(SIGMA)
                            Y_hat <- cbind(Y_hat, Y_temp)
						  } 
						}# -- End of "for 2" -- the last value of 'Y_temp' is the prediction T+N_pred_h
						Y_temp2[ii, ] <- Y_temp
					} # -- End of "for 1"
                #	Matrix of predictions               
                Y_pred[(((irep-N_burn)-1)*N_pred_thin+1):((irep-N_burn)*N_pred_thin), ] <- Y_temp2
                #	Predictive likelihood
                PL[(irep-N_burn), ] <- dmvnorm(Y1[T+N_pred_h, ], X_new_temp%*%ALPHA, SIGMA)
                if(PL[(irep-N_burn), ] == 0){
                    PL[(irep-N_burn), ] <- 1;
				  }
            }
        } # end forecasting
		#	=========Forecasting ends here

        #	=========IMPULSE RESPONSES:
        if(impulses == 1){
            biga <- matrix(0, M*p_lags, M*p_lags)
			j <- 1
            for(j in 1:(p_lags-1)){
                biga[(j*M+1):(M*(j+1)), (M*(j-1)+1):(j*M)] <- diag(M)
				}
            
            atemp <- ALPHA[2:dim(ALPHA)[1], ]
            atemp <- matrixcol(atemp)
            splace <- 0;
			ii <- 1
            for(ii in 1:p_lags){
                iii <- 1
				for(iii in 1:M){
                    biga[iii, ((ii-1)*M+1):(ii*M)] <- t(as.matrix(atemp[(splace+1):(splace+M),1]))
                    splace <- splace + M
					}
				}
            
            #	St dev matrix for structural VAR
            STCO <- chol(SIGMA)
            
            #	Now get impulse responses for 1 through nhor future periods
            impresp <- matrix(0, M, M*impulses_h)
            impresp[1:M,1:M] <- STCO
            bigai <- biga
			j <- 1
            for(j in 1:(impulses_h-1)){
                impresp[ ,(j*M+1):((j+1)*M)] <- bigj %*% bigai %*% t(bigj) %*% STCO
                bigai <- bigai%*%biga
				}
            
            #	Get the responses of all M variables to a shock imposed on
            #	the 'equatN'- th equation:
            equatN <- M													#	this assumes that the interest rate is sorted last in Y
            impf_m <- matrix(0, M, impulses_h)
            jj <- 0
			ij <- 1
            for(ij in 1:impulses_h){
                jj <- jj + equatN
                impf_m[ , ij] <- impresp[ ,jj]
				}
            imp[(irep-N_burn), , ] <- impf_m
		}
        #	=========IMPULSE RESPONSES ends here

    #	----- Save draws of the parameters
    alpha_draws[(irep-N_burn), ] <- t(alpha)
    ALPHA_draws[(irep-N_burn), , ] <- ALPHA
    SIGMA_draws[(irep-N_burn), , ] <- SIGMA
			
    }# -- End of "Store results"
		
}# -- End of "for 'irep'-- GIBBS iterations... "
#	=============================GIBBS SAMPLER ENDS HERE==================================
print('Time-processing elapsed:')
proc.time() - tic

#	Output
output <- list(alpha_draws,ALPHA_draws,SIGMA_draws)
return(output)
}
#
# -- End of "bvar.R" --