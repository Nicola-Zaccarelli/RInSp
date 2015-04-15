\name{MCp.RInSp}
\alias{MCp.RInSp}
\title{Monte Carlo resampling procedure}

\description{
The procedure replicates the two strategies of Monte Carlo resampling available in the following commands: \code{\link{Emc}}, \code{\link{PSicalc}}, \code{\link{WTdMC}}, and \code{\link{WTcMC}}.
}

\usage{
MCp.RInSp(dataset, pop.diet = "sum", replicates = 999)
}

\arguments{
  \item{dataset}{Object of class RInSp with data of type \dQuote{double} or \dQuote{integer}.}
  \item{pop.diet}{Option specifying the type of calculation for the population diet: \dQuote{sum} or \dQuote{average}. See \code{\link{pop.diet}} for details.}
  \item{replicates}{Number of resampled data matrices.}
}

\value{
A matrix of \dQuote{replicates + 1} datasets with the same number of individuals/resources of the input file. The first matrix is the original input data file.}
\details{
See Araujo et al (2011) for a discussion of the biological assumptions underlying this Monte Carlo resampling procedure. 

For counts of individual prey items representing approximately independent prey-capture decisions (i.e., integer data type), it is possible to run a Monte Carlo resampling simulation to test the null hypothesis that any observed diet variation arose from individuals sampling stochastically from a shared distribution. In the simulations, each individual is assigned a number of prey items equal to the number of items it was observed eating, and then prey items are randomly assigned to the individual's diet via multinomial sampling from the observed population resource distribution. Next, both \eqn{E} and \eqn{C_{ws}} are recalculated for the resulting simulated population.

For double data type, as mass of each food category in an individual's gut, the Monte Carlo resampling procedure first determines the number of prey used by each individual, \eqn{n_i}. The program then randomly reassigns each individual \eqn{n_i} prey items drawn from the population distribution of items. This is repeated a user-specified number of times. Hence, the resampled population yields a null model corresponding to a population composed of generalists that sample randomly from the population's diet, and have diet sizes equal to those of the observed data set.

The procedure can easily create a big three dimensional array. This object can create problems as it can consume an huge amount of memory. If needed use a package like \dQuote{ff} and a for loop to create an array able to hold the results.  
}
\references{
Araujo M.A., Layman C., and Bolnick D.I. 2011. The ecological causes of individual specialization. \emph{Ecology Letters} \bold{14}: 948-958.
}

\author{Dr. Nicola ZACCARELLI}

\seealso{
\code{\link{pop.diet}}, \code{\link{Emc}}, \code{\link{PSicalc}}, \code{\link{WTdMC}}, \code{\link{WTcMC}.}
}

\examples{
# Monte Carlo resampling for NODF
# EXAMPLE NOT RUN
# with stickleback data from Bolnick and Paull 2009
# data(Stickleback)
# Select a single spatial sampling site (site A)
# SiteA = import.RInSp(Stickleback, row.names = 1,
# info.cols = c(2:13), subset.rows = c("Site", "A"))
# nreplicates = 999
# Null.data = MCp.RInSp(SiteA, replicates = nreplicates)
# Build a vector holding real and resampled results
# the NODF values and warning messages from import.RInSp
# are printed, too.
# null.d.NODF = c()
# for (i in 1:(nreplicates + 1)) {
#   null.d.NODF = c(null.d.NODF, NODF(import.RInSp(Null.data[ , , i], print.messages=FALSE),
#  print.results=FALSE)$NODF)
# }
# hist(null.d.NODF)
# summary(null.d.NODF)
# Observed value of NODF from data 
# null.d.NODF[1]
# rm(list=ls(all=TRUE))
}
\keyword{univar}
\keyword{models}
