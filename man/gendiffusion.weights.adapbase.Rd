\name{gendiffusion.weights.adapbase}  
\alias{gendiffusion.weights.adapbase}
\title{Indices de Difusion de la Actividad Economica.}
\description{Esta funcion calcula cuatro indices de difusion de \emph{J} ciclos economicos desargrados y compuestos}
\usage{gendiffusion.weights.adapbase(niveles,crec_desagregado,crec_contribucion,
									 pesos,etiquetas)}
\arguments{
		\item{niveles}{Matriz de dimension \emph{(T X J)} con los niveles de los \emph{J} componentes de la actividad economica observados en \emph{T} periodos de tiempo.}
		\item{crec_desagregado}{Matriz de dimension \emph{(T-periodo) X J} con las variaciones relativas de los \emph{J} desagregados en \emph{(T-periodo)} observaciones temporales.}
		\item{crec_contribucion}{Matriz de dimension \emph{(T-periodo) X J} con las contribuciones de los \emph{J} desagregados en \emph{(T-periodo)} para las variaciones agregadas.}
		\item{pesos}{Matriz de dimension \emph{(A X J)} con los pesos de los \emph{J} desagregados para \emph{A} anios.}
		\item{etiquetas}{Etiquetas para las columnas de las entidades.}
		 }
\details{
         Inputs:
		 
		 \dQuote{niveles}:
		 
		 La primera columna del archivo debe contener el periodo (meses o trimestres). El archivo debe estar ordenado de manera ascendente respecto al periodo.
		
         La segunda columna del archivo \dQuote{niveles} debe contener el anio. Esta columna debe estar etiquetada con \dQuote{anio}.

		 \dQuote{pesos}:
		 
		 La primera columna del archivo debe contener el periodo (anio), este debe de estar ordenado de manera ascendente. 
		 		 
		 Las columnas restantes tanto de \dQuote{niveles} como de \dQuote{pesos} deben definirse con las entidades que se desee analizar como sigue:

		E1	.-	Entidad 1

		E2	.-	Entidad 2
		
			.
			
			.
			
			.
			
		EJ	.-	Entidad J
		}
\value{La funcion \dQuote{gendiffusion.weights.adapbase} genera como salidas una lista con los siguientes elementos:

			\item{agrgendiff}{Matriz de dimension \emph{(T-periodo) X 4} con los indices de difusion.}
			\item{diff0}{Convencional.}
			\item{diff1}{Ponderado (anio base).}
			\item{diff2}{Ponderado adaptativo.}
			\item{diff3}{Ponderado adaptativo y magnitud.}
Nota
			\item{periodo}{Hace referencia al periodo de diferenciacion que se aplica a los \emph{J} ciclos economicos.}
			\item{pesos}{Los calculos se realizan con diferentes pesos por anio hasta el anio mas reciente de acutalizacion de la base. A partir de ese anio, todos los claculos se realizan usando la base mas reciente.}
	   }

\references{
	\tabular{ll}{
		\tab - Guerrero-Escobar, S. y Martinez-Ovando, J. C. (2014) "Indices de Difusion de la Actividad Economica Estatal y Sectorial para Mexico." Sometido. \cr
		
		\tab - Guerrero-Escobar, S. y Martinez-Ovando, J. C. (2013) "A Note on the Diffusion of Business Cycles." Sometido.\cr
		}
	}
\author{
	\tabular{ll}{
				\tab Juan Carlos Martinez-Ovando <juan.martinez@banxico.org.mx>\cr

				\tab Santiago Guerrero-Escobar <sguerrero@banxico.org.mx>\cr
				}
		}
		
\examples{
	
	# Limpiamos la memoria 		
		rm(list = ls())

	# Cargamos la libreria

		library("gendifusion")

	# Definimos los inputs	
		
		data(datos)
			niveles <- datos
			data(pesos)
		etiquetas <- c("E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E12","E11","E13","E14",
					   "E15","E16","E17","E18","E19","E20","E21","E22","E23","E25","E24","E26",
					   "E27","E28","E29","E30","E31","E32")
		
		periodo <- as.integer(1)
		
#	Crecimiento desagregado
crec_desagregado <- growth.adapbase(niveles,periodo,etiquetas)
		
#-------------------------------------------------------------------------------------
#					crecimiento
#-------------------------------------------------------------------------------------

#	Crecimiento desagregado
crec_desagregado <- growth.adapbase(niveles,periodo,etiquetas)

#	Crecimiento agregado ponderado
crec_agr_salida <- weighted.growth.adapbase(niveles,pesos,periodo,etiquetas)
crec_agregado <- crec_agr_salida[[1]]
crec_contribucion <- crec_agr_salida[[2]]

#-------------------------------------------------------------------------------------
#					indices de difusion
#-------------------------------------------------------------------------------------

saldiff <- gendiffusion.weights.adapbase(niveles,crec_desagregado,crec_contribucion
										 ,pesos,etiquetas)
agrgendiff <- saldiff[[1]]
diff3_composicion <- saldiff[[2]]
		}			
		