repmat <- function(X,m,n){
  #
  #	R equivalent of repmat (matlab)
  #
  mx <- dim(X)[1]
  nx <- dim(X)[2]
  output <- matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
  #	Output
  return(output)
}
#
#	--  End of "repmat.R"