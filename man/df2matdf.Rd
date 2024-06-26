% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generic.R
\name{df2matdf}
\alias{df2matdf}
\title{Convert a response vector in a data frame to a matrix}
\usage{
df2matdf(x, formula)
}
\arguments{
\item{x}{a data frame}

\item{formula}{a formula}
}
\value{
A \code{data.frame}
}
\description{
Convert a response vector in a data frame to a matrix
}
\details{
This function identifies repeated combinations of explanatory variables
in a \code{mgcv} \code{formula} and then creates a \eqn{n \times m} 
\code{matrix} response variable in which each row corresponds to one of 
\eqn{n} unique explanatory variable combinations and each column to one of 
\eqn{m} replicates with the combination. Here \eqn{m} is the maximum number
of replicates for an explanatory variable combination; rows of the \code{matrix}
are padded with \code{NA}s at the end where there are fewer than \eqn{m}
replicates.
}
\references{
http://arma.sourceforge.net/docs.html#pinv
}
\seealso{
\link{match}, \link{unique}, \link{duplicated}
}
