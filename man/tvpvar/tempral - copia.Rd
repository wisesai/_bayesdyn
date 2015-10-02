\name{bvar}
\alias{bvar}
\title{Bayesian Vector Autoregression}
\description{Computes MCMC simulations from the posterior distribution of a Bayesian .}
\usage{tvpvar(ydata, yearlab, demean=0, N_training_prior, p_lags, N_mcmc, N_burn, N_print, prior_type=1)}
\arguments{
		\item{ydata}{ }
		\item{yearlab}{ }
		\item{demean}{ }
		\item{N_training_prior}{ }
		\item{p_lags}{ }
		\item{N_mcmc}{ }
		\item{N_burn}{ }
		\item{N_print}{ }
		\item{prior_type}{ }
        }        
\references{
		-	"Predictive Vector Autoregression under Instability and Informative Prior Distribution," J. C. Martinez-Ovando (2015).
		}
\author{
		Juan Carlos Martinez-Ovando <juan.martinez@banxico.org.mx>, <jc.martinez.ovando@gmail.com>
       }