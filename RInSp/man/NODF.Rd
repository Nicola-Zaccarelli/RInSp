\name{NODF}
\alias{NODF}
\title{
Nestedness metric based on overlap and decreasing fill}

\description{
The procedure calculates the \dQuote{nestedness metric based on overlap and decreasing fill} (NODF) for a binary matrix following Almeida-Neto et al. (2008).
}

\usage{
NODF(dataset, print.results= TRUE)
}

\arguments{
  \item{dataset}{Object of class RInSp with data of type \dQuote{double} or \dQuote{integer}.}
  \item{print.results}{Define if results for NODF should be printed. Default is TRUE}.
}

\value{
The result is a list of class \sQuote{RInSp} composed of:
  \item{NODF}{Value of the index of nestedness.}
  \item{Nrows}{Value of the index of nestedness for rows.}
  \item{Ncols}{Value of the index of nestedness for columns.}
  \item{R}{Binary matrix with individuals as rows and resources as columns. This matrix can be imported into the software PAJEK (\url{http://vlado.fmf.uni-lj.si/pub/networks/pajek/}) to draw a binary bipartite network of diet connectance between individuals (one set of nodes) and resources (a second set of nodes).}
  \item{NpR}{The degree of nestedness observed for rows.}
  \item{NpC}{The degree of nestedness observed for columns.}
}

\details{
Nestedness is a feature of binary matrices (also called presence/absence or incidence matrices). The notion of nestedness is particularly relevant for studies focusing the patterns of species occurrence among a set of locations (e.g., islands) and the patterns of interacting species within ecological networks. 

Almeida-Neto et al. (2008) propose a nestedness metric is based on two simple properties: decreasing fill (or DF) and paired overlap (or PO). Assuming that in a matrix with m rows and n columns, row i is located at an upper position from row j, and column k is located at a left position from column l. In addition, let MT be the marginal total (i.e. the sum of 1's) of any column or row. For any pair of rows/columns i and j, \eqn{DF_{ij}} will be equal to 100 if \eqn{MT_j} is lower than \eqn{MT_i}. Alternatively, \eqn{DF_{ij}} will be equal to 0 if \eqn{MT_j} is greater or equal to \eqn{MT_i}. For columns/rows, paired overlap (\eqn{PO_{kl}}) is simply the percentage of 1's in a given column/row l that are located at identical row/column positions to those in a column/row k. For any left-to-right column pair and, similarly, for any up-to-down row pair, there is a degree of paired nestedness (\eqn{N_{paired}}) as zero if \eqn{DF_{paired}} is zero, and PO if \eqn{DF_{paired}} is 100.

From the \eqn{n(n-1)/2}  and \eqn{m(m-1)/2} paired degrees of nestedness for n columns and m rows, we can calculate a measure of nestedness among all columns (\eqn{N_{col}}) and among all rows (\eqn{N_{row}}) by simply averaging all paired values of columns and rows.

Finally, the measure of nestedness for the whole matrix is given by:

\deqn{ NODF = \frac{\sum{N_{paired}}}{(\frac{n(n-1)}{2})+(\frac{m(m-1)}{2})} }

In the context of studies of individual specialization, one form of diet variation arises when individuals differ in their niche breadth, such that some individuals diet is a subset of other individuals' diet. This is revealed by a nestedness metric, which may be large (indicating nesting) or small (indicating clustering).
}

\references{
Almeida-Neto M., Guimaraes P., Guimaraes P.R., Loyola R.D., and Ulrich, W. 2008. A consistent metric for nestedness analysis in ecological systems: reconciling concept and measurement. \emph{Oikos} \bold{117}: 1227-1239.
}

\author{Dr. Nicola ZACCARELLI}

\seealso{For further indices see the R package \dQuote{vegan}. Please consider that the NODF implementation lacks the ordering and weighting options of \emph{nestednodf} available in \dQuote{vegan}.}

\examples{
# NODF example with stickleback data from Bolnick and Paull 2009
data(Stickleback)
# Select a single spatial sampling site (site D)
SiteD = import.RInSp(Stickleback, row.names = 1, info.cols = c(2:13),
subset.rows = c("Site", "D"))
Nesting = NODF(SiteD)
rm(list=ls(all=TRUE))
}

\keyword{univar}
\keyword{models}
