#--------------------------------------------------------
#	Uses of library: 	"dlm"
#--------------------------------------------------------

#
#	Author:		Juan Carlos Martinez-Ovando (juan.martinez@banxico.org.mx)
#

#             Read files with functions in R
rm(list=ls())

#	library("dlm") #	My package
library("MASS")
ruta = "K:/JCMO.Research/Coding/_r/_my_packages/_bayesdyn/bayesdyn/R/tvar"
source(paste(ruta,"/ar.R",sep=""))
source(paste(ruta,"/tvar.R",sep=""))
source(paste(ruta,"/tvar.pred.R",sep=""))
set.seed(12345)

#-------------------------------------------------------------
#		DATA
#-------------------------------------------------------------
rutadata = "K:/JCMO.Research/Coding/_r/_my_packages/_bayesdyn/bayesdyn/data"
datos <- read.csv(file = paste(rutadata,"/y.csv",sep=""), header = TRUE, sep = ",", quote = "\"",
         dec = ".", fill = TRUE, comment.char = "", row.names = 1)
datos <- as.data.frame(datos)
datos <- as.matrix(datos)
T <- nrow(datos)
x <- as.matrix(datos)
acf(x)

#	TVAR Model parameters
p <- 13
ar_fit <- ar(x,p)
m0 <- as.vector(ar_fit[[1]])
n0 <- 1
s0 <- 1
C0 <- 1*diag(p)
del <- c(0.994,0.995)

#	TVAR Model fit
tvar_fit <- tvar(x,p,del,m0,C0,s0,n0)
m <- tvar_fit[[1]]
C <- tvar_fit[[2]]
s <- tvar_fit[[3]]
n <- tvar_fit[[4]]
e <- tvar_fit[[5]]
ff <- tvar_fit[[6]]
fs <- tvar_fit[[7]]

x_ts <- ts(x, start = c(2001,1), frequency = 12)
e_ts <- ts(e, start = c(2001,1), frequency = 12)
ff_ts <- ts(ff, start = c(2001,1), frequency = 12) + mean(x)
fs_ts <- ts(fs, start = c(2001,1), frequency = 12) + mean(x)

plot(cbind(x_ts,ff_ts), plot.type = "single", lty = 1:3, col = 4:2, type = "b")   

#	TVAR Prediction
K <- 5
Msim <- 10000
tvar_pred <- tvar.pred(x,p,K,Msim,m,C,n,s)
tvar_pred_ts <- ts(tvar_pred, start = c(2001,1), frequency = 12)
plot(tvar_pred_ts, plot.type = "single", lty = 1:3, col = 4:2, type = "b")   

#	Fan chart
library('fanplot')
plot(NULL, type = "n", xlim = c(1, T+K), ylim = range(tvar_pred_ts), ylab = "Y", xlab = "Tiempo")
# add fan
fan(t(tvar_pred_ts))

#	--	END	--