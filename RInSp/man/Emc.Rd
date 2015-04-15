\name{Emc}
\alias{Emc}
\title{
Monte Carlo resampling procedure for index E and Cws
}
\description{
The procedure performs a Monte Carlo resampling under a null hypothesis to calculate the E measure of interindividual variation and the value of the \eqn{C_{ws}} measure of modularity following Araujo et al. (2008).
}
\usage{
Emc(dataset, popd.type = "sum", index = "saramaki", replicates=999)
}

\arguments{
  \item{dataset}{Object of class RInSp. Only data of type \dQuote{integer} can be used. See \code{\link{import.RInSp}} for details.}
  \item{popd.type}{Option specifying the way to calculate the population diet. Valid values are \dQuote{sum} or \dQuote{average}. See \code{\link{pop.diet}} for details.}
  \item{index}{The type of clustering coefficient to use. Valid values are Saramaki's and Barrat's index: \dQuote{saramaki} or \dQuote{barrat}. See \code{\link{Eindex}} for details.}
  \item{replicates}{Number of replicates for the Monte Carlo resampling, under the null hypothesis that any observed diet variation arose from individuals sampling stochastically from a shared distribution.}
}

\value{
The result is a list of class \sQuote{RInSp} composed of:
  \item{E}{Value of the index of individual specialization E.}
  \item{meannullE}{The mean value of the Monte Carlo resampling simulated data.}
  \item{Eadj}{The adjusted value of E.}
  \item{p.value}{It is the proportion of Monte Carlo resampling simulated values greater than the observed E.}
  \item{montecarlo}{It is a matrix where the first line collects the values of Omean, E, CW and Cws for the provided dataset, and the following lines hold the Monte Carlo resampled values.}
  \item{pop.diet}{Option specifying the way to calculate the population diet. Valid values are \dQuote{sum} or \dQuote{average}.}
  \item{ind.type}{The type of clustering coefficient to use. Valid values are Saramaki's and Barrat's index: \dQuote{saramaki} or \dQuote{barrat}.}
  \item{parameter}{an index used in \code{\link{sumMC.RInSp}} to point to the column of the matrix \sQuote{montecarlo} for plotting results. If zero no Monte Carlo resampling is performed.}
}

\details{
For counts of individual prey items representing approximately independent prey-capture decisions (i.e., integer data type), it is possible to run a Monte Carlo resampling simulation to test the null hypothesis that any observed diet variation arose from individuals sampling stochastically from a shared distribution. See Araujo et al. (2011) for a discussion of the biological assumptions underlying this Monte Carlo resampling procedure. In the simulations, each individual is assigned a number of prey items equal to the number of items it was observed eating, and then prey items are randomly assigned to the individual's diet via multinomial sampling from the observed population resource distribution. Next, both \eqn{E} and \eqn{C_{ws}} are recalculated for the resulting simulated population. The program generates a specified number of such null datasets, and the observed test statistic is compared to the distribution observed under the null hypothesis.

Note that the Monte Carlo method should not be applied to non-integer diet data (e.g., total mass of each prey species, proportion of stomach volume per prey species).

There are two options for calculating the population's diet proportions (\eqn{q_j}, the proportion of the resource j in the population's diet): \emph{sum} and \emph{average}.  When \emph{sum} is specified, one sums up all resource counts within a category across all individuals to get the population's use, then determine the proportion of each resource category in the population's repertoire:

\deqn{ q_j = \frac{\sum_i{n_{ij}}}{\sum_{ji}{n_{ij}}}}

The drawback of this approach is that individuals that eat large numbers of items, or larger total mass of items, will bias the population to look more like them.
The \emph{average} method (average proportion) circumvents this problem by first converting individual diets into proportions \eqn{p_{ik}}, then averaging these proportions for each resource k.

An adjusted value of E is calculated as:

\deqn{E_{adj} = \frac{E_{obs} - E_{null}}{1 - E_{null}}}

This rescales the value of E (i.e., the observed value \eqn{E_{obs}}) to range from a minimum of 0 when the observed value is equal to the mean Null value (i.e., \eqn{E_{null}}, up to a maximum of 1 when individual specialization is strongest). This allows one to compare datasets with different mean Null values.

Version 1.1 of the package fixes the case of highly specialised individuals in the calculation of \eqn{C_{ws}} by removing a bug in the C code taken from Dieta1.c affecting how the weight matrix and the number of links where calculated when highly specilised individuals are present.
}

\references{
Araujo M.S., Guimaraes Jr., P.R., Svanback, R., Pinheiro, A., Guimaraes P., dos Reis, S.F., and Bolnick, D.I. 2008. Network analysis reveals contrasting effects of intraspecific competition on individual vs. population diets. \emph{Ecology} \bold{89}: 1981-1993.

Araujo M.A., Layman C., and Bolnick, D.I. 2011. The ecological causes of individual specialization. \emph{Ecology Letters} \bold{14}: 948-958.

Schoener, T.W. 1968. The \emph{Anolis} lizards of Bimini: resources partitioning in a complex fauna. \emph{Ecology} \bold{49}: 704-726.
}

\author{Dr. Nicola ZACCARELLI.
Please note that the C code is partially based on the C source of Dieta.c by Araujo et al. (2008). The original source code Dieta1.c (identifier E089-115-A1; \url{http://esapubs.org/archive/ecol/E089/115/}) can be downloaded from the \dQuote{Ecological Archives} of the Ecological Society of America (\url{http://esapubs.org/archive/})}

\seealso{ Function \code{\link{Eindex}}, \code{\link{Null.Hp.RInSp}}, \code{\link{sumMC.RInSp}}.}

\examples{

# Example with stickleback data from Bolnick and Paull 2009
data(Stickleback)
# Select a single spatial sampling site (site B)
GutContents_SiteB = import.RInSp(Stickleback, row.names = 1,
info.cols = c(2:13), subset.rows = c("Site", "B"))
# Warning, the number of replicates is set low to speed up he example!
# Note, for real analyses we recommend to set replicates => 999
ResEmc = Emc(GutContents_SiteB, popd.type = "average", replicates = 249)

# This Eadj re-scales the observed E value to range from 0 when
# it is equal to the null expectation, up to a maximum of 1.
# This allows comparison of E values across samples with
# different null Es.
# For example, consider a case with no individual specialization
Nulldata = Null.Hp.RInSp(GutContents_SiteB, prop="average")
# Warning, the number of replicates is set low to speed up he example!
# Note, for real analyses we recommend to set replicates => 999
EmcNull = Emc(Nulldata, popd.type = "average", replicates = 249)
rm(list=ls(all=TRUE))
}

\keyword{univar}
\keyword{models}
