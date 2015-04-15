\name{WTdMC}
\alias{WTdMC}
\title{
Monte Carlo resampling of WIC/TNW for the discrete case
}
\description{
The program calculates the Total Niche Width (TNW), and breaks TNW down into its Between Individual Component (BIC) and Within Individual Component (WIC). It reports these three statistics and the proportion WIC/TNW. Niche widths of the population, within individuals, and among individuals are calculated using the Shannon-Weaver diversity index on count data, following Roughgarden's formulas (1974) and Bolnick et al. (2002).
A Monte Carlo resampling routine is used to calculate a null distribution for each of these statistics and generate a p-value testing the null hypothesis that all individuals sample equally from the population diet distribution. Note that the bootstrapping procedure is only biologically meaningful when applied to integer data representing counts of individual prey within individual predators' diets, and should not be applied to other forms of data such as prey mass or proportion of total volume.
}
\usage{
WTdMC(dataset, pop.diet = "sum", replicates=999, print.ris=TRUE)
}

\arguments{
  \item{dataset}{Object of class RInSp with data of type \dQuote{integer}.}
  \item{pop.diet}{Option specifying the way to calculate the population diet. Valid values are \dQuote{sum} or \dQuote{average}. See \code{\link{pop.diet}} for details.}
  \item{replicates}{Number of replicates for the Monte Carlo resampling.}
  \item{print.ris}{Printing results to screen. Default is TRUE.}
}

\value{
The result is a list of class \sQuote{RInSp} composed of:
\item{WonT}{The value of WIC/TNW for the dataset \dQuote{dataset}.}
\item{Zeros}{A character vector with the names of  individuals with Shannon-Weaver scores equal to zero. Or zero otherwise.}
\item{p.value}{The Monte Carlo resampling probability value for the calculated WIC/TNW ratio;}
\item{montecarlo}{A matrix with the replicate number \sQuote{Zeroes}, the value of \sQuote{WIC}, the value of \sQuote{BIC}, the value of \sQuote{TNW} and the derived value for \sQuote{WonT}. The first line holds the values for the provided dataset, and the following lines the results of the Monte Carlo resampling.}
\item{parameter}{An index used in \code{\link{sumMC.RInSp}} to point to the column of the matrix \sQuote{montecarlo} for plotting results.}
}

\details{
The Roughgarden's measures are the most widely reported indices of diet variation, due to historical precedent. They are closely analogous to an Analysis of Variance, in that one partitions the total variation in resource use (TNW) into within- and between-individual components (WIC and BIC). An advantage of this approach is that one can simultaneously evaluate the degree of diet variation (WIC/TNW), and shifts in individual and population niche widths. A drawback is that, for categorical diet data, the Shannon-Weaver measure of diversity can lead to biased estimates of diet variation (Bolnick et al. 2002). Despite its historical precedent, users are advised to avoid using WIC/TNW. This bias is not a concern for WIC/TNW measured on quantitative traits (prey size) rather than prey categories, but this is rarely done.

Roughgarden (1979, pg 510) proposed a measure to describe the niche width based on Shannon-Weaver information theory. The following formulae are equivalent to Roughgarden's formulation, although we use different notation. Let N be a matrix of diet data, where elements \eqn{n_{ij}} represent the number (or mass) of diet items in individual i's diet that fall in category j. This raw data matrix is then transformed into a proportion matrix P, with elements \eqn{p_{ij}} describing the proportion of the jth resource category in individual (or phenotype) i's diet. 

Using these proportions or the raw numerical counts one then calculates population diet proportions (\eqn{q_{k}}, see \code{\link{pop.diet}} for details). These data can then be used to calculate WIC, BIC, and TNW (Roughgarden 1979; Bolnick et al. 2002):

\deqn{WIC = \sum_i{p_i} (- \sum_k{p_{ik} ln(p_{ik})})}

\deqn{BIC = \sum_k{p_i ln(p_i)} - \sum_k{q_k (- \sum_i{t_{ik}ln(t_{ik})})}}

\deqn{TNW = - \sum_k{q_k ln(q_k)}}

where: \eqn{p_i} is the proportion of all resources used by individual i (\eqn{\frac{\sum_k{n_{ik}}}{\sum_{ik}{n_{ik}}}}); \eqn{q_k} is the proportion of the kth resource category in the population's niche, and \eqn{t_{ik}} is the proportion of the population's total use of resource k that was used by individual i (i.e., \eqn{\frac{n_{ik}}{\sum_i{n_{ik}}}}).

It is important to remember that TNW = BIC + WIC, and one can easily calculate the proportion WIC, which we refer to as WIC/TNW.

Be aware that monophagous individuals (100\% of diet from a single resource) have excessively large leverage drawing WIC towards zero, possibly leading to downward-biased WIC/TNW exaggerating individual specialization. The procedure prints a warning message with the name of monophagous individuals in the dataset.

The Monte Carlo resampling procedure first determines the number of prey used by each individual, \eqn{n_i}. It then determines the frequencies with which the population as a whole used the various resources (see \code{\link{pop.diet}} for details). Subsequently, the procedure randomly reassigns each individual \eqn{n_i} prey drawn from the population frequency distribution, and then recalculates a null measure of diet variation. This is repeated a user-specified number of times. Hence, the resampled population yields a null model corresponding to a population composed of generalists that sample randomly from the population's diet, and have diet sizes equal to those of the observed data set.
The original statistics are merged into the Monte Carlo set.
}

\references{
Bolnick, D.I., L.H. Yang, J.A. Fordyce, J.M. Davis, and Svanback, R. 2002. Measuring individual-level resource specialization. \emph{Ecology} \bold{83}: 2936-2941.

Roughgarden, J. 1979. Theory of population genetics and evolutionary ecology: an introduction. Macmillan, New York, NY, USA. 
}
\author{Dr. Nicola ZACCARELLI. Please note that the C code is partially based on the C++ source of IndSpec from Bolnick et al. (2002).

The original source code of IndSpec1.0.exe is available from the \dQuote{Ecological Archives} of the Ecological Society of America (\url{http://esapubs.org/archive/}): identifier E083-056-S1; \url{http://esapubs.org/archive/ecol/E083/056/}.}

\seealso{ For the continuous case see \code{\link{WTcMC}} and for a summary statistics the function \code{\link{sumMC.RInSp}}.}

\examples{
data(Stickleback)
# Select a single spatial sampling site (site A)
GutContents_SiteA = import.RInSp(Stickleback, row.names = 1,
info.cols = c(2:13), subset.rows = c("Site", "A"))
# Warning, the number of replicates is set low to speed up he example!
# Note, for real analyses we recommend to set replicates => 999
WT = WTdMC(GutContents_SiteA, replicates = 99)
sumMC.RInSp(WT)
rm(list=ls(all=TRUE))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{univar}
\keyword{models}
