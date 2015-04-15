\name{overlap}
\alias{overlap}
\title{Mean pairwise overlap}
\description{
This procedure calculates the pairwise diet overlap between all individuals in a sample.
}
\usage{
overlap(dataset)
}
\arguments{
  \item{dataset }{Object of class RInSp with data.}
  }
\value{
Return an matrix of three columns:
  \item{meanoverlap}{The mean pairwise overlap for all pairs of individuals, excluding the diagonal.}
  \item{meanindividualoverlap}{The mean pairwise overlap for individual i and all others, excluding the diagonal.}
  \item{meandissimilarity}{The mean pairwise dissimilarity overlap for all pairs of individuals, excluding the diagonal.}
  \item{overlapmatrix}{The probability associated to the likelihood value.}
  \item{parameter=0}{Internal option for plotting results.}
}
\details{
For a sample size of N individuals, this will produce an N * N matrix, each cell \eqn{o_{ik}} representing the diet overlap between individual i and individual k. The diagonal is all ones because an individual has a 100\% overlap with itself. Pairwise overlap is calculated by determining the proportional similarity between each individual.

The script reports the average pairwise diet overlap (not including diagonal values), the average pairwise diet dissimilarity (i.e., one minus the average pairwise diet overlap) and the complete diet overlap matrix. The overlap matrix can be useful for testing whether diet (dis)similarity is a function of various other metrics of between-individual difference or similarity (e.g., Bolnick and Paull 2009).}

\references{
Bolnick, D.I., and Paull, J.S., 2009. Morphological and dietary differences between individuals are weakly but positively correlated within a population of threespine Stickleback. \emph{Evolutionary Ecology Research} \bold{11}: 1217-1233.

Bolnick, D.I., L.H. Yang, J.A. Fordyce, J.M. Davis, and Svanback, R. 2002. Measuring individual-level resource specialization. \emph{Ecology} \bold{83}: 2936-2941.
}

\author{Dr. Nicola ZACCARELLI}
\examples{
 data(Stickleback)
# Select a single spatial sampling site (site D)
GutContents_SiteD = import.RInSp(Stickleback, row.names = 1,
info.cols = c(2:13), subset.rows = c("Site", "D"))
similarity = overlap(GutContents_SiteD)

\dontrun{
# Generate a plot similar to Bolnick and Paull 2009, testing whether
# morphological differences between individuals are positively related
# to diet differences (or conversely are negatively related to diet similarity)
GutContents = import.RInSp(Stickleback, row.names = 1,info.cols = c(2:13))
similarity = overlap(GutContents)
diet.similarity = as.dist(similarity$overlapmatrix)
size.dissimilarity = dist(GutContents$info[,3])
plot(diet.similarity ~ size.dissimilarity, xlab = "Size dissimilarity",
ylab = "Diet similarity", pch = 16, cex = 0.4)
model = lm(diet.similarity ~ size.dissimilarity)
abline(model, lwd = 3, col = "red")

# The ade4 package is needed to complete the analysis
library(ade4)
mantel.rtest(diet.similarity, size.dissimilarity)
}
}
\keyword{function}
\keyword{utilities}
