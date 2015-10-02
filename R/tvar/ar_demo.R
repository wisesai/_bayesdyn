#--------------------------------------------------------
#	Uses of library: 	"dlm"
#--------------------------------------------------------

#             Read files with functions in R
rm(list=ls())

#	library("dlm") o funciones
#library("matrixcalc")
ruta = "L:/JCMO.Research/Coding/_r/_my_packages/_bayesdyn/bayesdyn/R/tvar"
source(paste(ruta,"/ar.R",sep=""))
set.seed(12345)

#-------------------------------------------------------------
#		DATA
#-------------------------------------------------------------
rutadata = "I:/JCMO.Research/_banxico/proyectos/_frutas"
datos <- read.csv(file = paste(rutadata,"/indice_frutas_verduras_00_13.csv",sep=""), header = TRUE, sep = ",", quote = "\"",
         dec = ".", fill = TRUE, comment.char = "", row.names = 1)
datos <- as.data.frame(datos)
T <- nrow(datos)
x <- as.matrix(datos[14:T,"diflog"])

#	Model parameters and AR(p) fit
p <- 8
ar_fit <- ar(x,p)

#	Decomposition
m <- as.vector(ar_fit[[1]])

#	--	END	--