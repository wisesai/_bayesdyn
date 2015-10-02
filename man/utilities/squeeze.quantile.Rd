\name{squeeze.quantile}
\alias{squeeze.quantile}
\title{Squeeze quantiles of three dimensional arrays}
\description{Computes quantiles of a three-dimensional array (dim JxKxL) with respect to the first dimension and the probability vector 'probs'. Output is an array of dimension KxLxlength(probs).}
\usage{squeeze.quantile(Array,probs)}
\arguments{
		\item{Array}{ }
		\item{probs}{ }
        }        
\references{
		-	"Predictive Vector Autoregression under Instability and Informative Prior Distribution," J. C. Martinez-Ovando (2015).
		}
\author{
		Juan Carlos Martinez-Ovando <juan.martinez@banxico.org.mx>, <jc.martinez.ovando@gmail.com>
       }