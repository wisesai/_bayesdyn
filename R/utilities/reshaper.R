reshaper <- function(A, type){
#
#	This function reshapes the matrix A of dim MxN into 
#	vectors of two types:
#		1	-	matrix of dime MN x 1 (by columns) 
#		2	-	matrix of dime 1 x NM (by rows) 
#	
#	Author:
# 	Juan~Carlos Martínez-Ovando
#	Email:		juan.martinez@banxico.org.mx
#				jc.martinez.ovando@gmail.com
#

M <- dim(A)[1]		#	num. of rows
N <- dim(A)[2]		#	num. of columns

if(type == 'col'){
	Arep <- matrix(0,M*N,1)
	count <- 0
	k <- 1
	for(k in 1:N){
		Arep[(count*M + 1):(count*M + M), 1] <- A[, k]
		count <- count + 1
	  }#	-- End of "for 1"
  }else if(type == 'row'){
	Arep <- matrix(0,1,M*N)
	count <- 0
	j <- 1
	for(j in 1:M){
		Arep[1, (count*N + 1):(count*N + N)] <- A[j, ]
		count <- count + 1
	  }#	-- End of "for 1"
  }

#	Output
return(Arep)
}
#
#	--	End of "reshaper.R"