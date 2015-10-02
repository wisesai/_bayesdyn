\name{dlmbma}
\alias{dlmbma}
\title{Bayesian Dynamic Model Averaging (Filtering and Smoothing)}
\description{This function computes Bayesian model avraging, forward-filtering and backward smoothing of a general homoscedastic regression dynamic linear model among several model confiburations.}
\usage{dlmbma(x, y, models.which, lambda=0.99, gammaa=0.99, eps=.001/nrow(models.which), modelthreshold=21)}
\arguments{
		\item{x}{ }
		\item{y}{ }
		\item{models.which}{ }
		\item{lambda}{ }
		\item{gammaa}{ }
		\item{eps}{ }
		\item{modelthreshold}{ }
        }        
\references{
    - "A Quantitative Assessment of the Influence of Crime on the Economic Growth Dynamics in Mexico in the Short-Term," Martinez-Ovando, J. C. (2014).
		}
\author{
		Juan Carlos Martinez-Ovando <juan.martinez@banxico.org.mx>, <jc.martinez.ovando@gmail.com>
       }