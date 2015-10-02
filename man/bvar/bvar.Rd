\name{bvar}
\alias{bvar}
\title{Bayesian Vector Autoregression}
\description{Computes MCMC simulations from the posterior distribution of a Bayesian VAR model. The function also computes simulations from the predictive distribution of the VAR model in two fashion (direct or iterated), and impulse-response to shocks.}
\usage{bvar(Yraw, constant, p_lags, forecasting, forecast_method, N_pred_h, impulses, impulses_h, prior, N_mcmc, N_burn, N_print)}
\arguments{
		\item{Yraw}{ }
		\item{constant}{ }
		\item{p_lags}{ }
		\item{forecasting}{ }
		\item{forecast_method}{ }
		\item{N_pred_h}{ }
		\item{impulses}{ }
		\item{impulses_h}{ }
		\item{prior}{ }
		\item{N_mcmc}{ }
		\item{N_burn}{ }
		\item{N_print}{ }
        }        
\references{
		-	"Predictive Vector Autoregression under Instability and Informative Prior Distribution," J. C. Martinez-Ovando (2015).
		}
\author{
		Juan Carlos Martinez-Ovando <juan.martinez@banxico.org.mx>, <jc.martinez.ovando@gmail.com>
       }