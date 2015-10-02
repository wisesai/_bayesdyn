\name{gibbs_ss}
\alias{gibbs_ss}
\title{Gibbs sampler for state-space models}
\description{This function implements the Gibbs sampler for state space models based on the algorithm proposed by Carter and Kohn (1994).}
\usage{gibbs_ss(y,Z,Ht,Qt,m,p,t,B0,V0)}
\arguments{
		\item{y}{ }
		\item{Z}{ }
		\item{Ht}{ }
		\item{Qt}{ }
		\item{m}{ }
		\item{p}{ }
		\item{t}{ }
		\item{B0}{ }
		\item{V0}{ }
        }        
\references{
		-	"Predictive Vector Autoregression under Instability and Informative Prior Distribution," J. C. Martinez-Ovando (2015).
		}
\author{
		Juan Carlos Martinez-Ovando <juan.martinez@banxico.org.mx>, <jc.martinez.ovando@gmail.com>
       }