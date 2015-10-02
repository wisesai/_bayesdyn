reshapes <- function(Arep, type, M, N){
#
#	reshapes.R
#	This function reshapes an array A of dim (MN)x1 (or 1x(MN)) into 
#	matrices of two types:
#		1	-	matrix of dime M x N (by columns in 'Arep') 
#		2	-	matrix of dime M x N (by rows in 'Arep') 
#	
#	Note:	This function complements "reshaper.R"
#
#	Author: 	Juan~Carlos Martínez-Ovando
#	Email:		juan.martinez@banxico.org.mx
#				jc.martinez.ovando@gmail.com
#

if(length(Arep) != M*N ){print('Error in reshapes.R')}

if(type == 'col'){
	A <- matrix(0,M,N)
	count <- 0
	k <- 1
	for(k in 1:N){
		A[, k] <- Arep[(count*M + 1):(count*M + M), 1] 
		count <- count + 1
		}#	-- End of "for 1"
}else if(type == 'row'){
	A <- matrix(0,M,N)
	count <- 0
	j <- 1
	for(j in 1:M){
		A[j, ] <- Arep[(count*N + 1):(count*N + N), 1]
		count <- count + 1
		}#	-- End of "for 1"
}

return(A)
}
#
#	--	End of "reshapes.R"