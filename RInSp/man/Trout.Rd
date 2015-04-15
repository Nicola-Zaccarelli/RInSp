\name{Trout}
\alias{Trout}
\docType{data}
\title{Example of continuous data from fish prey lengths of brown trout}
\description{
The example reproduces data published in Kahilainen & Lehtonen (2001) on prey lengths from a sample of 59 individuals of stocked (S) and native (N) brown trout (\emph{Salmo trutta}) in the subarctic Lake Muddusjarvi in northern Finland.

Column names refer to:

- FishID: a unique identifier for each individual fish within the geographic area;

- Type: stocked (S) and native (N) brown trout.

All other columns are lengths (in mm) of fish prey found in stomachs.
}
\usage{data(Trout)}

\references{
Kahilainen, K. and Lehtonen, H. 2001. Resource use of native and stocked brown trout \emph{Salmo trutta} L., in a subarctic lake. \emph{Fisheries Management and Ecology} \bold{8}: 83-94.
}

\examples{
# Summary of total length of preys in stomach content
# by geographic region
data(Trout)
# Mean prey length by individual
# First change zeros to NA
troutTMP = Trout
troutTMP[ troutTMP == 0] = NA
individuals = by(troutTMP[, 3:7], troutTMP[, 2], rowMeans, na.rm=TRUE)
TroutN = subset(Trout, Type == "N", PL1:PL5)
TroutN = TroutN[TroutN > 0]
TroutS = subset(Trout, Type == "S", PL1:PL5)
TroutS = TroutN[TroutS > 0]
boxplot(c(TroutS, TroutN) ~ c(rep("N", length(TroutN)), rep("S", length(TroutS))))
rm(list=ls(all=TRUE))

}
\keyword{dataset}
