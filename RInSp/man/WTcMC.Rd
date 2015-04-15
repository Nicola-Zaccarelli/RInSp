\name{WTcMC}
\alias{WTcMC}
\title{
Monte Carlo resampling of WIC/TNW for the continuous case
}
\description{
The program calculates the Total Niche Width (TNW), and breaks TNW down into its Between Individual Component (BIC) and Within Individual Component (WIC). It reports these three statistics and the proportion WIC/TNW. Niche widths of the population, within individuals, and among individuals are calculated following Roughgarden's formulas (1974).

A Monte Carlo resampling routine is used to calculate a null distribution for each of these statistics and generate a p-value testing the null hypothesis that all individuals sample equally from the population diet distribution.
}

\usage{
WTcMC(dataset, replicates = 999, weight = "equal", print.ris=TRUE)
}

\arguments{
  \item{dataset}{Object of class RInSp with data of type \dQuote{double}.}
  \item{replicates}{Number of replicates for the Monte Carlo resampling.}
  \item{weight}{The weighting scheme for individuals. Default value is \dQuote{equal}. Alternative value is \dQuote{N_items}.}
  \item{print.ris}{Printing results to screen. Default is TRUE.}
}

\value{
The result is a list of class \dQuote{RInSp}composed of:
  \item{WonT}{The value of WIC/TNW for the dataset \dQuote{dataset}.}
  \item{p.value}{The Monte Carlo resampling probability value for the calculated WIC/TNW ratio.}
  \item{montecarlo}{is a matrix with the replicate number \sQuote{REP}, the value of \sQuote{WIC}, the value of \sQuote{BIC}, the value of \sQuote{TNW} and the derived value for \sQuote{WonT}.  The first line holds the values for the provided dataset, and the following lines the results of the Monte Carlo resampling.}
  \item{weight}{The weighting scheme for individuals used.}
  \item{parameter}{An index used in \code{\link{sumMC.RInSp}} to point to the column of the matrix \sQuote{montecarlo} for plotting results.}
}

\details{
Roughgarden (1972) suggested that the total niche width of a population (TNW) can be broken down into two components: the variation in resource use within individuals (within-individual component, WIC), and the variance between individuals (between-individual component, BIC) so that TNW = WIC + BIC. Originally  developed  for  within-  and  between-phenotype (or between-sex) variation, it can be applied at the individual level. Assuming that variation in niche parameters can be expressed along a single continuous dimension x (e.g., prey size), let X be a matrix of diet data where each element \eqn{x_{ij}} is the size (or other measure) of the jth prey item in individual i's diet. Then,

\deqn{TNW = Var(x_{ij})}

\deqn{WIC = E(Var(x_j | i))}

\deqn{BIC = Var(E(x_j | i))}

The relative degree of individual specialization can be measured as the proportion of TNW explained by within-individual variation, WIC/TNW. As this value approaches 1, all individuals utilize the full range of the population's niche, whereas smaller values indicate decreasing inter-individual overlap and hence higher individual specialization.

Roughgarden's WIC/TNW is limited to continuous diet data. To carry this approach over to discrete data such as the frequency of alternate prey taxa in the diet, Roughgarden (1979:510) proposed a measure that uses the Shannon-Weaver index as a proxy for variance (See \code{\link{WTdMC}} for details).

In case of a marked presence of individuals with one or few prey items, the relation TNW = BIC + WIC may not hold. In this case there is the option to weight the number of items per individuals following two alternatives:

- \dQuote{equal}: weighting by the number of items in each individual's diet, so those with more data contribute more to estimating parameters; \cr
- \dQuote{N_items}: weighting each individual equally regardless of diet items number; requires weighting each diet item by the inverse of the number of items in the individual's diet for calculating TNW.

The Monte Carlo resampling procedure first determines the number of prey used by each individual, \eqn{n_i}. Subsequently, the procedure randomly reassigns each individual \eqn{n_i} prey items drawn from the population distribution of items. This is repeated a user-specified number of times. Hence, the resampled population yields a null model corresponding to a population composed of generalists that sample randomly from the population's diet, and have diet sizes equal to those of the observed data set.
The original statistics are merged into the Monte Carlo set.
}

\references{
Bolnick, D.I., L.H. Yang, J.A. Fordyce, J.M. Davis, and Svanback, R. 2002. Measuring individual-level resource specialization. \emph{Ecology} \bold{83}: 2936-2941.

Ingram,T., Stutz, W.E., Bolnick, D.I. 2011. Does intraspecific size variation in a predator affects its diet diversity and top-down control of prey? \emph{PLoS ONE} \bold{6}: e20782.

Roughgarden, J. 1974. Niche width: biogeographic patterns among \emph{Anolis lizard} populations. \emph{American Naturalist} \bold{108}: 429-411. 
}

\author{Dr. Nicola ZACCARELLI. Please note that the C code is partially based on the C++ source of IndSpec from Bolnick et al. (2002). For the \dQuote{weight} option Travis Ingram is acknowledged (Ingram et al. 2011).

The original source code of IndSpec1.0 is available from the \dQuote{Ecological Archives} of the Ecological Society of America (\url{http://esapubs.org/archive/}): identifier E083-056-S1; \url{http://esapubs.org/archive/ecol/E083/056/}.}

\seealso{ Function \code{\link{sumMC.RInSp}}.}

\examples{
# Example of analysis with the trout data from Finland
data(Trout)
TroutRIS= import.RInSp(Trout, col.header=TRUE, row.names=1, info.cols=2, data.type="double")
# Warning, the number of replicates is set low to speed up he example!
# Note, for real analyses we recommend to set replicates => 999
TroutRIS_MC = WTcMC(TroutRIS, replicates = 99)
sumMC.RInSp(TroutRIS_MC)
rm(list=ls(all=TRUE))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{univar}
\keyword{models}
