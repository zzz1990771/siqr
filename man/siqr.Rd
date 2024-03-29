% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/siqr.R
\name{siqr}
\alias{siqr}
\title{Main estimation function of single index quantile regression model.
a two step method.}
\usage{
siqr(y, X, tau = 0.5, beta.initial = NULL, h = NULL, maxiter = 30, tol = 1e-08)
}
\arguments{
\item{y}{response vector;}

\item{X}{covariate matrix;}

\item{tau}{left-tail probability (quantile index), scalar}

\item{beta.initial}{starting value of beta, the single index coefficients}

\item{h}{user-defined bandwidth}

\item{maxiter}{max iteration number}

\item{tol}{toleration for convergence}
}
\value{
a siqr object, which includes:
        beta - the fitted single index coefficients with unit norm and first component being non negative
        flag.conv  - whether the iterations converge
}
\description{
Main estimation function of single index quantile regression model.
a two step method.
}
\examples{
#generate data
set.seed(2021)
data <- generate.data(50)
X <- data$X
y0<- data$Y

#initials
beta0 <- NULL
#quantile
tau = 0.75
siqr.result <- siqr(y0,X,beta.initial = beta0, tau=tau)
summary(siqr.result)

}
