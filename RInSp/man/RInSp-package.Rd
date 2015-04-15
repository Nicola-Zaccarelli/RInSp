\name{RInSp-package}
\alias{RInSp-package}
\alias{RInSp}
\docType{package}
\title{
R Individual Specialization
}
\description{
\dQuote{R Individual Specialization} (RInSp) is an R package for analysing resource use data to calculate ecological measures of within-population niche variation based on data on individuals' resource use (e.g., stomach contents, repeated observations of prey choice or habitat use). The package implements several ecological niche metrics to measure individual and population niche width in order to measure individual specialization. This version of RInSp is able to compute Araujo's E, measures of clustering and pairwise similarity among individuals, IS, Petraitis' W, and Roughgarden's WIC/TNW, and an example of a measure of nestedness (NODF). For some of these indices Monte Carlo re-sampling procedures for testing significance are provided.

To cite package check the command \code{CITATION} or use:

Zaccarelli, N., Mancinelli, G. and Bolnick, D.I. (2013) RInSp: an R package for the analysis of individual specialisation in resource use. Methods in Ecology and Evolution, 4(11): 1018-1023, doi: 10.1111/2041-210X.12079

This package can be used on Mac, Windows and UNIX platforms, and on both i386 and x64 architectures. It is based on IndSpec1.0 developed by Prof. D. Bolnick (Bolnick et al. 2002), and on Dieta1 developed by Dr. M.A. Araujo (Araujo et al. 2008).

The original source code of IndSpec1.0 is available from the \dQuote{Ecological Archives} of the Ecological Society of America (\url{http://esapubs.org/archive/}): identifier E083-056-S1; \url{http://esapubs.org/archive/ecol/E083/056/}.
The original source code Dieta1 is available from the \dQuote{Ecological Archives} of the Ecological Society of America (\url{http://esapubs.org/archive/}): identifier E089-115-A1; \url{http://esapubs.org/archive/ecol/E089/115/}.

Version 1.1 of the package fixes three issues:
\enumerate{
\item in \code{\link{Eindex}} the case of highly specialised individuals is now handled correctly;
\item in \code{\link{like.Wi}} the probability for each value of Petraitis' W is properly handled;
\item in \code{\link{Emc}} the C code taken from Dieta1.c has been fixed for an error in the calculation of the weight matrix when highly specilised individuals are present.}

The Authors thank for the first two issues Miss Georgia Karoline Kosmala (\email{geh_karol@hotmail.com}) and  Miss Ambika Kamath (\email{geh_karol@hotmail.com}), and for starting the discussion leading to the identification of the bug in Dieta1.c Dr. Marcio S. Araujo (\email{msaraujo@rc.unesp.br}).

}

\details{
\tabular{ll}{
Package: \tab RInSp\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2015-02-10\cr
License: \tab GPL (>= 2)\cr
}
}
\author{
Dr. Nicola ZACCARELLI \email{nicola.zaccarelli@gmail.com}, \cr
Dr. Giorgio MANCINELLI \email{giorgio.mancinelli@unisalento.it}, \cr
and Prof. Dan BOLNICK \email{danbolnick@mail.utexas.edu}}

\references{
Almeida-Neto M., Guimaraes P., Guimaraes P.R., Loyola R.D., and Ulrich W. 2008. A consistent metric for nestedness analysis in ecological systems: reconciling concept and measurement. \emph{Oikos} \bold{117}: 1227-1239.

Araujo M.S., Guimaraes Jr., P.R., Svanback, R., Pinheiro, A., Guimaraes P., dos Reis, S.F., and Bolnick, D.I. 2008. Network analysis reveals contrasting effects of intraspecific competition on individual vs. population diets. \emph{Ecology} \bold{89}: 1981-1993.

Araujo M.A., Layman C., and Bolnick D.I. 2011. The ecological causes of individual specialization. \emph{Ecology Letters} \bold{14}: 948-958.

Aversen. J.N. 1969. Jackknifing U-statistics. \emph{Annals of Mathematical Statistics} \bold{40}: 2076-2100.

Barrat A., Barthelemy M., Pastor-Satorras R., and Vespignani A. 2004. The architecture of complex weighted networks. \emph{Proceedings of the National Academy of Sciences} \bold{101}: 3747-3752

Bolnick, D.I., L.H. Yang, J.A. Fordyce, J.M. Davis, and Svanback, R. 2002. Measuring individual-level resource specialization. \emph{Ecology} \bold{83}: 2936-2941.

Bolnick, D.I., and Paull, J.S. 2009. Morphological and dietary differences between individuals are weakly but positively correlated within a population of threespine Stickleback. \emph{Evolutionary Ecology Research} \bold{11}: 1217-1233.

Onnela J.P., Saramaki J., Kertesz J., and Kaski K. 2005. Intensity and coherence of motifs in weighted complex networks. \emph{Physics Review E} \bold{71}: 065103.

Roughgarden, J. 1974. Niche width: biogeographic patterns among \emph{Anolis lizard} populations. \emph{American Naturalist} \bold{108}: 429-411.

Saramaki J., Kivela M., Onnela J.P., Kaski K., and Kertesz J. 2007. Generalizations of the clustering coefficient to weighted complex networks. \emph{Physics Review E} \bold{75}: 027105.

Schoener, T.W. 1968. The \emph{Anolis} lizards of Bimini: resources partitioning in a complex fauna. \emph{Ecology} \bold{49}: 704-726.
}

