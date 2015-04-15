\name{import.RInSp}
\alias{import.RInSp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Import data into a RInSp object
}
\description{
The procedure reads and checks data to create an object of class RInSp.
}
\usage{
import.RInSp(filename, col.header = FALSE, row.names = 0, info.cols = 0,
subset.column = 0, subset.rows = NA, data.type = "integer", print.messages=TRUE)

}
\arguments{
  \item{filename}{Name of the file or dataframe to be read.}
  \item{col.header}{Logical value to indicate if a header row is there. Default is FALSE.}
  \item{row.names}{Numeric value for column of rows' name. Default is zero for no names.}
  \item{info.cols}{A vector collecting columns numbers for additional information. Defaults is to have no additional information.}
  \item{subset.column}{A vector of columns' indices to be used as a subset. Default is to have no columns subsetting.}
  \item{subset.rows}{A string vector where the first element points to the column name to be used for rows subsetting, and following elements with criteria. Default is to have no rows subsetting.}
  \item{data.type}{Data type among "integer" or "double". From the data a proportion matrix will be produced. Default is to use integers/counts values.}
  \item{print.messages}{Prints messages concerning the number of rows and columns eventually deleted after subsetting because composed of all zeros. Default is TRUE.}
}

\details{
Three different types of data can be used. Integer/count values and decimal/real values are stored in the \dQuote{resources} section of the output list and used to derive all meaningful information. While for proportions the \dQuote{resources} section will be empty. Valid key words are: \dQuote{integer}, \dQuote{double}, and \dQuote{proportion}. Use zero for empty cells.

The procedure will check for the presence of zero sum columns/rows for the selected dataset. A warning is printed in case of column/row deletion. It must be considered that the degree of checking on the subsetting is low.
}

\value{
Return an list of class RInSp composed by:
  \item{resources}{A matrix of the resources data.}
  \item{proportions}{A matrix of proportions of each resources in its row. This matrix can be imported into the software PAJEK (\url{http://vlado.fmf.uni-lj.si/pub/networks/pajek/}) to draw a weighted bipartite network connecting individuals to the various prey categories.}
  \item{data.type}{Data type used.}
  \item{col.names}{Name of the different columns of the resource data.}
  \item{ind.names}{Name of the individual.}
  \item{info}{A data frame containing additional information for the resource data.}
  \item{num.prey}{Number of resources/prey (i.e., columns) in the dataset after zero sum checking.}
  \item{num.individuals}{Number of individuals/sites (i.e., rows) in the dataset after zero sum checking.}
  \item{num.zero.prey}{Number of resources/prey (i.e., columns) in the dataset without zero sum checking.}
  \item{num.ind.zero}{Number of individuals/sites (i.e., rows) in the dataset without zero sum checking.}
}

\author{Dr. Nicola ZACCARELLI}
\examples{

data(Stickleback)
# Import data
GutContents = import.RInSp(Stickleback,  row.names = 1, info.cols = c(2:13))
GutContents
# Select a single spatial sampling site (site A)
GutContents_SiteA = import.RInSp(Stickleback, row.names = 1, info.cols = c(2:13),
subset.rows = c("Site", "A"))

# Select a subset of prey types
GutContents_subset = import.RInSp(Stickleback, row.names = 1, info.cols = c(2:13),
subset.column = c(13:28, 45))

# Lump prey types into functional groups then import data
# define new columns representing lumped prey categories
attach(Stickleback)
Copepods = Calanoid + Cyclopoid + Harpacticoid
Diptera = Diptera.Pupae + Chironomid.larvae + Ceratopogonid.larvae + Tipulid.larvae +
Tipulidae.Adult + Diptera.Adult + Diptera.Larvae + Ceratopogonid.Adult
InsectLarvae = Ephemeroptera + Trichoptera.larvae + Ephemeroptera.pupae +
Zygoptera.larvae + Plecoptera.larvae
Cladocera = Bosmina + Polyphemus + Holopedium + Daphnia + Chydorus
names(Stickleback)
GutContents_lumped = import.RInSp(Stickleback, row.names = 1, info.cols = c(2:13),
subset.column = c(18,24,28,31,43:46))
rm(list=ls(all=TRUE))
}
\keyword{function}
\keyword{utilities}
