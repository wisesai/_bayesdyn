\name{reshaper}
\alias{reshaper}
\title{Reshape a matrix array into a vector array}
\description{This function reshapes the matrix A (dim MxN) into vectors of two types: 1	-- matrix of dim MN x 1 (by columns), 2	-- matrix of dim 1 x NM (by rows).}
\usage{reshaper(A, type)}
\arguments{
		\item{A}{ }
		\item{type}{Type of transformation: 'col' (by columns) or 'row' (by rows)}
        }        
\references{
		-	"Predictive Vector Autoregression under Instability and Informative Prior Distribution," J. C. Martinez-Ovando (2015).
		}
\author{
		Juan Carlos Martinez-Ovando <juan.martinez@banxico.org.mx>, <jc.martinez.ovando@gmail.com>
       }