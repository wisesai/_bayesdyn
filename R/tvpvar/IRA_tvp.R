        if(istore == 1){
            
        #    
        #    IRA_tvp	% in Matlab
        #
 		    #	Impulse response analysis. Note that Htsd contains the
            #	structural error cov matrix
            #	Set up things in VAR(1) format as in Lutkepohl (2005) page 51
            biga <- matrix(0, M*p_lags,M*p_lags)
			j <- 1
            for(j in 1:(p_lags-1)){
                biga[(j*M+1):(M*(j+1)),(M*(j-1)+1):(j*M)] <- diag(M)
				}

            i <- 1
			for(i in 1:t){										#	Get impulses recurssively for each time period
                bbtemp <- as.matrix(Btdraw[(M+1):K,i])			#	get the draw of B(t) at time i=1,...,T  (exclude intercept)
                splace <- 0
				ii <- 1
                for(ii in 1:p_lags){
					iii <- 1
                    for(iii in 1:M){
                        biga[iii,((ii-1)*M+1):(ii*M)] <- t(as.matrix(bbtemp[(splace+1):(splace+M),1]))
                        splace <- splace + M
						}# -- End of "for 3"
					}# -- End of "for 2"

                #	------------Identification code:                
                #	St dev matrix for structural VAR
                Hsd <- Htsd[((i-1)*M+1):(i*M),1:M]				#	First shock is the Cholesky of the VAR covariance
                diagonal <- diag(diag(Hsd))
                Hsd <- solve(diagonal) %*% Hsd					#	Unit initial shock					

                #	Now get impulse responses for 1 through nhor future periods
                impresp <- matrix(0, M, M*nhor)
                impresp[1:M, 1:M] <- Hsd						#	First shock is the Cholesky of the VAR covariance
                bigai <- biga
				j <- 1
                for(j in 1:(nhor-1)){
                    impresp[ ,(j*M+1):((j+1)*M)] <- bigj %*% bigai %*% t(bigj) %*% Hsd
                    bigai <- bigai %*% biga
					}

                #	Only for specified periods
				yearlab <- as.matrix(yearlab)
                if(yearlab[i,1] == 1975.00){					#	1975:Q1
                    impf_m <- matrix(0, M, nhor)
                    jj <- 0
					ij <- 1
                    for(ij in 1:nhor){
                        jj <- jj + M							#	restrict to the M-th equation, the interest rate
                        impf_m[ ,ij] <- as.matrix(impresp[ ,jj])
						}
                    imp75[irep-N_burn, , ] <- impf_m				#	store draws of responses
					}
                if(yearlab[i,1] == 1981.50){					#	1981:Q3
                    impf_m <- matrix(0, M, nhor)
                    jj <- 0
					ij <- 1
                    for(ij in 1:nhor){
                        jj <- jj + M							#	restrict to the M-th equation, the interest rate
                        impf_m[ ,ij] <- impresp[ ,jj]
						}
                    imp81[(irep-N_burn), , ] <- impf_m			#	store draws of responses
					}
                if(yearlab[i,1] == 1996.00){					#	1996:Q1
                    impf_m <- matrix(0, M, nhor)
                    jj <- 0
					ij <- 1
                    for(ij in 1:nhor){
                        jj <- jj + M							#	restrict to the M-th equation, the interest rate
                        impf_m[ ,ij] <- impresp[ ,jj]
						}
                    imp96[(irep-N_burn), , ] <- impf_m			#	store draws of responses
					}
						
				}# -- End of "for 1" "geting impulses for each time period"
				
        } # -- End of "impulse response calculation section"  