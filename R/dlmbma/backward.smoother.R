backward.smoother <- function(m0,C0,n0,s0,Ft,Gt,yt){
#	
#	This function implements the backward-smoother procedure for
#	individual dinamic linear models
#
#	INPUTS:
#
#
#	OUTPUTS:
#
#
#	Author: 	Juan~Carlos Martínez-Ovando
#	Email:		juan.martinez@banxico.org.mx
#				JC.Martinez.Ovando@gmail.com
#
#	Reference:	Martínez-Ovando, J.~C. (2014) "The Influence of Crime on the Economic Growth 
#					Dynamics in Mexico in the Short-Term," Banco de México, Mimeo.
#

#	One-step backward
at <- Gt %*% m0
Rt <- t(Gt) %*% C0 %*% Gt #+ Wt
ft <- Ft %*% at
Qt <- drop(Ft %*% Rt %*% Ft + s0)
#	Updating...
et <- drop(yt-ft)
At <- (Rt %*% Ft) / Qt
nt <- n0 + 1
st <- drop(s0 + (s0/nt)*((et^2)/Qt - 1))
mt <- at + At * et
colnames(mt) <- colnames(m0)
rownames(mt) <- rownames(m0)
Ct <- (st/drop(s0)) * (Rt - (At%*%t(At)) * Qt)
colnames(Ct) <- rownames(mt)
rownames(Ct) <- rownames(mt)

# Output
list(ft=ft,
	 Qt=Qt,
	 mt=mt, 
	 Ct=Ct, 
	 nt=nt, 
	 st=st)
}
#
#	END of "forward.filter.R"