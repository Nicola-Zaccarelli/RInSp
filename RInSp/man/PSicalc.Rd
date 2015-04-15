\name{PSicalc}
\alias{PSicalc}
\title{Calculate the proportional similarity index}
\description{
A measure of individual specialization proposed by Bolnick et al. (2003), based on the average pairwise overlap of the niche distribution of individuals and the population.
}

\usage{
PSicalc(dataset, pop.diet = "sum", exclude = FALSE, replicates=999)
}
\arguments{
  \item{dataset}{Object of class RInSp with data.}
  \item{pop.diet}{Option specifying the way to calculate the population diet. Valid values are \dQuote{sum} or \dQuote{average}.}
  \item{exclude}{A logical value to exclude or not the individual in variance estimates. Defalut is FALSE, i.e. keep all individuals.}
  \item{replicates}{Number of replicates for Monte Carlo resampling.}
}

\value{
The function returns a list composed by:
 \item{PSi}{A column vector with the individual proportional similarity index for each individual i.}
 \item{IS}{The value of the prevalence of individual specialization.}
 \item{PSi.montecarlo}{Monte Carlo resampling simulation values for \eqn{PS_i}}
 \item{Var.montecarlo}{Monte Carlo resampling simulation values for \eqn{PS_i}'s variance}
 \item{VarPSi}{A column vector of the estimated variance for each individual.}
 \item{population.diet}{A vector describing the population diet. A matrix is returned when \dQuote{exclude} is TRUE.}
 \item{IS.pvalue}{Probability value for the observed IS against the Monte Carlo resampling simulation. Consider to use more than 1000 simulations.}
 \item{montecarlo}{Monte Carlo resampling simulation values for \eqn{PS_i}}
 \item{num.individuals}{The number of individuals of the data set.}
 \item{parmeter}{The value points to the column of the \dQuote{montecarlo} matrix for plotting results.} 
 }
\details{
Both Schoener (1968) and Feinsinger et al. (1981) advocated using Czekanowski's proportional similarity index (PS) for interspecific niche studies. Adapted to individual-level analyses, the diet overlap between an individual and the population is: 

\deqn{ PS_i = 1 - 0.5 \sum_j{|p_{ij} - q_j|} }

where \eqn{p_{ij}} is the frequency of category j in the individual i's diet, and \eqn{q_j}  is the frequency of category j in the population as a whole.
The prevalence of individual specialization (IS) in the population is then measured by the average of individuals' PS values:

\deqn{IS = \frac{\sum_i(PS_i)}{N}} 

Monte Carlo resampling for hypothesis testing is implemented for the case when all individuals are kept.

The user has the option of calculating each individual's \eqn{PS_i} and \eqn{Var(PS_i)} using \eqn{q_j} values calculated by either including or excluding individual i.}
\references{
Bolnick, D.I., L.H. Yang, J.A. Fordyce, J.M. Davis, and Svanback, R. 2002. Measuring individual-level resource specialization. \emph{Ecology} \bold{83}: 2936-2941.

Feinsinger, P., Spears, E. E., and Poole, R. W. 1981. A simple measure of niche breadth. \emph{Ecology}  \bold{62}: 27-32.

Schoener, T.W. 1968. The \emph{Anolis} lizards of Bimini: resource partitioning in a complex fauna. \emph{Ecology} \bold{49}: 704-726.
}

\author{Dr. Nicola ZACCARELLI}

\examples{
# PSicalc example using stickleback data from Bolnick and Paull 2009
data(Stickleback)
# Select a single spatial sampling site (site A)
GutContents_SiteA = import.RInSp(Stickleback, row.names = 1,
info.cols = c(2:13), subset.rows = c("Site", "A"))
# Warning, the number of replicates is set low to speed up he example!
# Note, for real analyses we recommend to set replicates => 999
PSi = PSicalc(GutContents_SiteA, exclude = FALSE, replicates = 99)
sumMC.RInSp(PSi)
rm(list=ls(all=TRUE))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{function}
\keyword{utilities}
