\name{Stickleback}
\alias{Stickleback}
\docType{data}
\title{
Example of data from threespine Stickleback gut contents
}
\description{
Example of a raw data matrix to be used with the \code{\link{import.RInSp}} procedure. Data are from Bolnick and Paull (2009), studying individual specialization in threespine Stickleback (\emph{Gasterosteus aculeatus}) from 5 sites within one lake in British Columbia. Data columns include:

- FishID: a unique identifier for each individual fish;

- Site: A-E, five locations around the lake where fish were collected;

- Mass: fish mass in grams;

- Standard.length, Body.width, gape.width: fish dimensions in centimetres;

- raker.number: gill raker number;

- raker.length: gill raker length (mm);

- Sex: sex;

- MicrogN, DeltaAir, MicrogC, DeltaPDB: stable isotope data for C and N stable isotope ratios.

All other columns are counts of prey categories in Stickleback stomach contents.

}
\usage{data(Stickleback)}

\references{
Bolnick, D.I., and Paull, J.S. 2009. Morphological and dietary differences between individuals are weakly but positively correlated within a population of threespine Stickleback. \emph{Evolutionary Ecology Research} \bold{11}: 1217-1233.
}

\examples{
data(Stickleback)
summary(Stickleback)
rm(list=ls(all=TRUE))
}
\keyword{dataset}
