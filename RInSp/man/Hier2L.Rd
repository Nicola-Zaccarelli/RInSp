\name{Hier2L}
\alias{Hier2L}
\title{
Variance partition and WIC/TNW
}
\description{
Variance partition based on a single multi-level factor and comparison of the relative contribution of the \dQuote{\emph{Between Group Component}} (BGC) versus the ratio WIC/TNW. For a more general discussion of variance partition the reader is pointed to Viola et al. (2012).
}

\usage{
Hier2L(dataset, factor = 1, weight.type = "N_items")
}

\arguments{
  \item{dataset}{Object of class RInSp with data of type \dQuote{double}.}
  \item{factor}{Numeric value indicating the column in the \emph{info} part of the RInSp input dataset containing the factor used for build subsets.}
  \item{weight.type}{The weighting option for unbalanced prey numbers. See details in \code{\link{WTcMC}}.}
}

\value{
The result is a matrix of 5 rows by n+1 columns, when n is the number of levels of the factor used to create groups in the dataset. The first four rows are the values for WIC, BIC, TNC and WIC/TNW for all possible levels and the whole dataset. Last row is the value of the \dQuote{\emph{Between Group Component}} (BGC) for the n-levels factor.
}

\details{
The Beteen Group Component is calculated as the (population) variance of the mean values of the n-levels groups.

The level of checking performed by the procedure is low. Pay attention to input parameters.
}

\references{

Violle C., Enquist B.J., McGill B.J., Jiang L., Albert C.H., Hulshof C., Jung V. and Messier J. 2012. The return of the variance: intraspecific variability in community ecology. \emph{Trends in Ecology & Evolution} \bold{27}: 244252.
}

\author{Dr. Nicola ZACCARELLI and Dr. Giorgio MANCINELLI}

\seealso{Function \code{\link{WTcMC}}.}

\examples{
# Example of variance decomposition
data(Trout)
TroutRIS= import.RInSp(Trout, col.header=TRUE, row.names=1, info.cols=2, data.type="double")
decomp = Hier2L(TroutRIS, factor=1)
rm(list=ls(all=TRUE))
}

\keyword{univar}
\keyword{models}
