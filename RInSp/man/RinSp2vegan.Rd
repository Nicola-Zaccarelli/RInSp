\name{RInSp2vegan}
\alias{RInSp2vegan}
\title{
Vegan support procedure}

\description{
The procedure converts an RInSp class object into a dataframe to be used as input for the \dQuote{Vegan} package to calculate nestedness indices.
}
\usage{
RInSp2vegan(dataset)
}

\arguments{
  \item{dataset}{Object of class RInSp with data of type \dQuote{double}, \dQuote{integer} or \dQuote{proportions}.}
}

\value{
The result is a binary dataframe with rows and columns names from the input object.
}

\author{Dr. Nicola ZACCARELLI}

\seealso{For further indices see the R package \dQuote{vegan}.}

\examples{
data(Stickleback)
# Select a single spatial sampling site (site D)
GutContents = import.RInSp(Stickleback, row.names = 1, info.cols = c(2:13))
risul = RInSp2vegan(GutContents)
rm(list=ls(all=TRUE))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{univar}
\keyword{models}
