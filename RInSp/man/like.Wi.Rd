\name{Likelihood and Wi}
\alias{like.Wi}
\title{Likelihood measure of niche breadth}
\description{
The procedure calculates the estimation of likelihood measures of niche breadth and overlap described in Petraitis (1979).}
\usage{
like.Wi(dataset)
}
\arguments{
  \item{dataset }{Object of class RInSp with data.}
  }

\value{
Return a list of class RInSp with:
  \item{MeanWi}{the mean population value of Wi;}
  \item{ResCat}{the number of resource categories;}
  \item{ind.vals}{A matrix with three columns: \dQuote{Likelihood} with value of the likelihood index for the individual \emph{i}; \dQuote{p-value} for the the probability associated to the likelihood value; \dQuote{Wi} with the value of the Petraitis' W index for the individual \emph{i}.}
}
\details{
The function returns the likelihood of the observed diet (\eqn{\lambda_i}) the associated probability , and the value of the Petraitirs' W.
The likelihood of the observed diet of individual \emph{i} is:

\deqn{\lambda_i = \prod_j (\frac{q_j}{p_{ij}})^{n_{ij}}}

where \eqn{q_j} is the population proportion of the resource \emph{j}, \eqn{p_{ij}} is the proportion of the resource \emph{j} in the diet of the individual \emph{i}, and \eqn{n_{ij}} is the number of items for the individual \emph{i} and the resource \emph{j}.

This can be used to calculate a p-value to test the significance of the diet specialization, as \eqn{-2ln(\lambda)} is distributed as a chi-square with (r-1) degrees of freedom, where r is the number of resource categories.

The generalised likelihood ratio test rejects the null hypothesis for a unilateral alternative hypotesis using significance level \eqn{\alpha} if:

\deqn{-2ln(\lambda) > \chi^2_{(r-1)}}

Petraitis' W is computed following:
\deqn{W_i = \lambda_i^{(1/D_i)}}
where \eqn{D_i} is the number of diet items recorded in the diet of individual \emph{i}. This index is a measure of niche width relative to a specified distribution. For a complete generalist individual, \eqn{W_i = 1}, and the value decreases with greater specialization.
}

\references{
Petraitis, P. S. 1979. Likelihood measures of niche breadth and overlap. \emph{Ecology} \bold{60}(4): 703-710. 

Bolnick, D.I., L.H. Yang, J.A. Fordyce, J.M. Davis, and Svanback, R. 2002. Measuring individual-level resource specialization. \emph{Ecology} \bold{83}: 2936-2941.
}
\author{Dr. Nicola ZACCARELLI}

\examples{
# Likelihood and Wi example with stickleback data
# from Bolnick and Paull 2009
data(Stickleback)
# Select a single spatial sampling site (site D)
SiteD = import.RInSp(Stickleback, row.names = 1,
info.cols = c(2:13), subset.rows = c("Site", "D"))
Wi = like.Wi(SiteD)
rm(list=ls(all=TRUE))
}
\keyword{function}
\keyword{utilities}
