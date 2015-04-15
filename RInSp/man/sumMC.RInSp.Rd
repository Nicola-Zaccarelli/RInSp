\name{sumMC.RInSp}
\alias{sumMC.RInSp}
\title{
Summary function for Monte Carlo resampling results
}
\description{
The function  provides summary statistics on the Monte Carlo resampling simulated index distributions. It plots an histogram with the position of the actual value of the index for the original dataset (blue line) and the 2.5\% and 97.5\% percentile of the simulated distribution (red dashed lines).

It works for objects of class \sQuote{RInSp}.
}
\usage{
sumMC.RInSp(dataset)
}

\arguments{
  \item{dataset}{Dataset of class RInSp.}
}

\author{Dr. Nicola ZACCARELLI}

\seealso{ See also \code{\link{WTcMC}}, \code{\link{WTdMC}}, or \code{\link{Emc}}.}

\examples{
#  sumMC.RInSp example using stickleback data from Bolnick and Paull 2009
data(Stickleback)
# Select a single spatial sampling site (site A)
GutContents_SiteA = import.RInSp(Stickleback, row.names = 1,
info.cols = c(2:13), subset.rows = c("Site", "A"))
EmcRes = Emc(GutContents_SiteA, replicates=499)
sumMC.RInSp(EmcRes)
rm(list=ls(all=TRUE))
}
\keyword{function}
\keyword{utilities}
