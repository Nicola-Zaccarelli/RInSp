\name{pop.diet}
\alias{pop.diet}
\title{
Calculate population diet
}
\description{
This function calculates the population diet using two different approaches: the \dQuote{sum} or the \dQuote{average} method.
}

\usage{
pop.diet(dataset, prop = "sum")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{Object of class RInSp with data.}
  \item{prop}{Option specifying the type of calculation. Two types: \dQuote{sum} or \dQuote{average}.}
}

\value{
Returns a list of information:
 \item{popdiet}{Vector of the specified population diet.}
 \item{popdtype}{The type of population diet.}
 \item{richness}{The richness of resources (for now it has sense only for integer data type).}
 \item{D}{The Levins' D index for resources.}
}
\details{
The \dQuote{sum} approach is the most straightforward. One sums up all resource counts within a category across all individuals to get the population's use, then determine the proportion of each resource category in the population's repertoire. The proportion \eqn{q_j} of the resource j in the population's diet is:

\deqn{q_j = \frac{\sum_i{n_{ij}}}{\sum_{ji}{n_{ij}}}}

The drawback of this approach is that individuals that eat large numbers of items, or larger total mass of items, will bias the population to look more like them.

The \dQuote{average} method (average proportion) circumvents this problem by first converting individual diets into proportions \eqn{p_{ik}}, then averaging these proportions for each resource k.

Along with the population's diet the procedure calculates the Levins' D index (Levins 1968) of diversity as:

\deqn{ D = 1 - \frac{1}{\sum{q_{j}^2}}}
}
\references{
Bolnick, D.I., L.H. Yang, J.A. Fordyce, J.M. Davis, and Svanback, R. 2002. Measuring individual-level resource specialization. \emph{Ecology} \bold{83}: 2936-2941.

Levins R. 1968. \emph{Evolution in Changing Environments: Some Theoretical Explorations}. Princeton Univ. Press.
}
\author{Dr. Nicola ZACCARELLI}

\seealso{ Function \code{\link{Null.Hp.RInSp}}.}

\examples{
data(Stickleback)
# Select a single spatial sampling site (site D)
GutContents_SiteD = import.RInSp(Stickleback, row.names = 1, info.cols = c(2:13),
subset.rows = c("Site", "D"))
PopDiet = pop.diet(GutContents_SiteD, prop = "average")
rm(list=ls(all=TRUE))
}

\keyword{function}
\keyword{utilities}
