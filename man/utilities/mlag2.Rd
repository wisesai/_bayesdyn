\name{mlag2}
\alias{mlag2}
\title{Generates a matrix with lagged values}
\description{Generates a matrix of n lags from a matrix (or vector) containing a set of vectors -- for use in VAR model routines.}
\usage{mlag2(Y, p)}
\arguments{
		\item{Y}{A matrix (or vector), nobs x nvar}
		\item{p}{Number of contiguous lags for each vector in x}
        }        
\references{
		-	"Predictive Vector Autoregression under Instability and Informative Prior Distribution," J. C. Martinez-Ovando (2015).
		}
\author{
		Juan Carlos Martinez-Ovando <juan.martinez@banxico.org.mx>, <jc.martinez.ovando@gmail.com>
       }