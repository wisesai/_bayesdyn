\name{reshapes}
\alias{reshapes}
\title{Reshape a vector array into a matrix array}
\description{This function reshapes an array A of dim (MN)x1 (or 1x(MN)) into matrices of two types: 1	-- matrix of dim M x N (by columns in 'Arep'),  2	-- matrix of dim M x N (by rows in 'Arep'). Note: This function complements "reshaper.R"}
\usage{reshapes(A, type,  M, N)}
\arguments{
		\item{A}{ }
		\item{type}{Type of transformation: 'col' (by columns) or 'row' (by rows)}
		\item{M}{ }
		\item{N}{ }
        }        
\references{
		-	"Predictive Vector Autoregression under Instability and an Informative Prior Distribution," J. C. Martinez-Ovando (2015).
		}
\author{
		Juan Carlos Martinez-Ovando <juan.martinez@banxico.org.mx>, <jc.martinez.ovando@gmail.com>
       }