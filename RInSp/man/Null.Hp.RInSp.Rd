\name{Null.Hp.RInSp}
\alias{Null.Hp.RInSp}
\title{
Generate a data set based on a specific null hypothesis}

\description{
The procedure let you build an integer/count dataset of class RInSp:

- by fixing the number of items/resources per individual and the population diet;

- by providing a specific dataset of class RInSp and and option of \dQuote{sum} or \dQuote{average} for the population diet.

It builds a table with a specified number of individuals and for each, a resources vector is created using a multinomial distribution with probability equal to the specified population diet and number of items.

As the procedure checks for zero sum columns/rows, it can be time consuming.
}

\usage{
Null.Hp.RInSp(dataset, prop = "sum")
}

\arguments{
  \item{dataset}{Object of class RInSp with data of type \dQuote{integer} or a numeric integer vector specifying the number of items for each individual in the simulated data set.}
  \item{prop}{Option specifying the type of calculation for population diet. Two types are valid when \dQuote{dataset} is an object of class RInSp: \dQuote{sum} or \dQuote{average}. For details see \code{\link{pop.diet}}. Or a numeric vector when the counts of items for individuals is provided.}
}

\value{
An object of class RInSp holding the null model data set.
}

\details{
The procedure can take a long time as it assure that the final data set has no columns/rows sum equal to zero.

It can be used under two scenarios:

- when a specific data set of integer/count data is of interest, it can be used as input to generate a null data set with the same number of individuals/resources and total number of diet items, but the probability of assignment can be derived by the population diet (providing the \dQuote{sum} or \dQuote{average}) or a diet specified by the user;

- to generate a null data set by specifying the total number of items for each individual and the probability of the resource type in the diet using a diet provided.

When the user provides the number of items for individuals, the procedure checks for negative values and it converts the numbers into integers.

When the diet is specified by the user, the procedure will convert all values into positive numbers and recalculate the proportions of each entry, so to ensure to have a vector of probabilities (i.e., sum equal to one, values between zero and one).

A progress bar shows the progress in generating random samples till a maximum of 400 new random data sets.

The procedure can easily create a big three dimensional array. This object can create problems as it can consume an huge amount of memory. If needed use a package like \dQuote{ff} and a for loop to create an array able to hold the results.
}

\author{Dr. Nicola ZACCARELLI and Dr. I. Bolnick}

\seealso{ Function \code{\link{pop.diet}}.}

\examples{
# Example with stickleback data from Bolnick and Paull 2009
data(Stickleback)
# Null model using a specified data set
SiteA = import.RInSp(Stickleback, row.names = 1,
info.cols = c(2:13), subset.rows = c("Site", "A"))
Eresult = Eindex(SiteA, index = "saramaki", jackknife = FALSE)
Nulldata = Null.Hp.RInSp(SiteA, prop="average")
Eresult = Eindex(Nulldata, index = "saramaki", jackknife = FALSE)

# Null model with diet and per individual resource
# Generate a number-of-items-per-individual vector for 30 individuals
ind.resources = floor(runif(30, 15, 120))
# Generate a population diet with 6 items
diet = c(10, 30, 15, 60, 100)
Nulldata2 = Null.Hp.RInSp(ind.resources, prop= diet)
Eresult2 = Eindex(Nulldata2, index = "saramaki", jackknife = FALSE)
rm(list=ls(all=TRUE))
}

\keyword{univar}
\keyword{models}
