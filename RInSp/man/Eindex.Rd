\name{Eindex}
\alias{Eindex}
\title{
Calculate index E and Cws
}
\description{
The procedure calculates the E measure of interindividual variation, its variance and the value of the \eqn{C_{ws}} measure of modularity after Araujo et al. (2008).
}

\usage{
Eindex(dataset, index = "saramaki", jackknife = FALSE)
}

\arguments{
  \item{dataset}{Object of class RInSp with data of type \dQuote{double}, \dQuote{integer} or \dQuote{proportions}.}
  \item{index}{The type of clustering coefficient to use. Valid values are Saramaki's and Barrat's index: \dQuote{saramaki} or \dQuote{barrat}.}
  \item{jackknife}{Specify if a jackknife estimate of the index variance is required. Default is FALSE.}
}

\value{
The result is a list of class \sQuote{RInSp} composed of:
  \item{Omean}{Mean value of the measure of the network overall degree of pairwise overlap.}
  \item{E}{Value of the index of individual specialization E.}
  \item{PS}{Matrix of the measure of niche pairwise overlap between i and j, adapted from Schoener (1968).}
  \item{PSbinary}{Binary matrix derived by applying the threshold of Omean to the PS matrix. This matrix can be imported into the software PAJEK (\url{http://vlado.fmf.uni-lj.si/pub/networks/pajek/}) to draw binary networks of diet similarity among individuals (e.g., see Araujo et al. 2008).}
  \item{Ejack}{Values of the measure of interindividual variation used for the Jackknife estimates of Var(E).}
  \item{VarE}{Jackknife estimates of the variance of E.}
  \item{CW}{is the network weighted clustering coefficient.}
  \item{CwS}{Vector of the individuals weighted clustering coefficients.}
  \item{Cw}{Value of the measure of modularity.}
  \item{index}{The type of clustering coefficient used.}
  \item{Ki}{A vector with the degree of the nodes of the network.}
}

\details{

The index E has been proposed by Araujo et al. (2008) as a measure of individual specialization where, in absence of interindividual niche variation, its value is zero. The index will increase towards one with the increase of interindividual variation.

A jackknife estimation of the variance of E can be derived using the formalism of U-statistics (Aversen, 1969). For a complete description and a formal demonstration of the jackknife estimation of the variance of the E index the reader is refered to Araujo et al. (2008).

A measure of the relative degree of clustering in a network to test for modularity in the niche overlap network is \eqn{C_{ws}}. In a totally random network (i.e., a network consisting of individuals that sample randomly from the population niche), \eqn{C_{ws}} is approximally 0, indicating no modularity. If individuals form discrete groups specialized on distinct sets of resources, \eqn{C_{ws} > 0}, and the network is modular. If \eqn{C_{ws} < 0} the network degree of clustering is actually lower than what would be expected solely on the overall network density of connections, indicating that diet variation takes place at the level of the individual, as opposed to discrete groups. 

The relative degree of clustering \eqn{C_{ws}} is obtained as the mean values over all nodes the niche overlap network of the individual node weighted clustering coefficients (\eqn{C_{w_i}}).

The clustering can be measured using two different type of weighted clustering coefficients. In general the degree of unweighted clustering around a vertex \emph{i} in a network is quantified by evaluating the number of triangles in which the vertex participates normalized by the maximum possible number of such triangles. Hence we have zero if none of the neighbours of a vertex are connected, otherwise we have one. By extending the above line of reasoning, the weighted clustering coefficient should also take into account how much weight is present in the neighbourhood of the vertex, compared to some limiting case. This can be done in several ways.
Barrat et al. (2004) were the  first to  propose a weighted version of the clustering coefficient of the form:

\deqn{ C_{w_i} = \frac{1}{s_i(k_i -1)} \sum_{j,h}{\frac{(w_{ij}+w_{ih})}{2}a_{ij}a_{ih}a_{jh}} }

where \eqn{s_i} is the sum of the weights (\eqn{w_i}) of all the edges between node i and the nodes to which it is connected; \eqn{k_i} is the number of edges between node i and its neighbours; \eqn{w_{ij}} is the weight of the edge between thw two nodes i and j; \eqn{a} is 1 if an edge is present between each pair \emph{ij}, \emph{ih}, and \emph{jh} respectively, and zero otherwise. The summation, therefore, quantifies the weights of all edges between node i and its neighbours that are also neighbours to each other.

Samaraki et al. (2007) and Onnela et al. (2005) proposed a version of a clustering index of the form:

\deqn{ C_{w_i} = \frac{1}{k_i(k_i -1)} \sum_{j,h}{(w_{ij}w_{ih}w_{jh})^\frac{1}{3}} }

where \eqn{k_i} is the number of edges between individual i and its neighbours; \eqn{w_{ij}} is the weight of the edge between individual i and j obtained by dividing the actual weight by the maximum of all weights. The summation, therefore, quantifies the weights of all edges between individual i and its neighbours that are also neighbours to each other.

The default value of the procedure has been set for analogy with the Dieta1.c code provided by Araujo et al. (2008). The original source code Dieta1.c is available from the \dQuote{Ecological Archives} of the Ecological Society of America (\url{http://esapubs.org/archive/}): identifier E089-115-A1; \url{http://esapubs.org/archive/ecol/E089/115/}.

Version 1.1 of the package fixes the case of highly specialised individuals in the calculation of \eqn{C_{ws}}.
}

\references{
Araujo M.S., Guimaraes Jr., P.R., Svanback, R., Pinheiro, A., Guimaraes P., dos Reis, S.F., and Bolnick, D.I. 2008. Network analysis reveals contrasting effects of intraspecific competition on individual vs. population diets. \emph{Ecology} \bold{89}: 1981-1993.

Aversen J.N. 1969. Jackknifing U-statistics. \emph{Annals of Mathematical Statistics} \bold{40}: 2076-2100.

Barrat A., Barthelemy M., Pastor-Satorras R., and Vespignani, A. 2004. The architecture of complex weighted networks. \emph{Proceedings of the National Academy of Sciences} \bold{101}: 3747-3752

Onnela J.P., Saramaki J., Kertesz J., and Kaski, K. 2005. Intensity and coherence of motifs in weighted complex networks. \emph{Physics Review E} \bold{71}: 065103.

Saramaki J., Kivela M., Onnela J.P., Kaski K., and Kertesz, J. 2007. Generalizations of the clustering coefficient to weighted complex networks. \emph{Physics Review E} \bold{75}: 027105.

Schoener, T.W. 1968. The \emph{Anolis} lizards of Bimini: resources partitioning in a complex fauna. \emph{Ecology} \bold{49}: 704-726.
}

\author{Dr. Nicola ZACCARELLI}

\seealso{ Function \code{\link{Emc}}.}

\examples{
# Eindex example with data from Bolnick and Paull (2009)
data(Stickleback)
# Select a single spatial sampling site (site B)
GutContents_SiteB = import.RInSp(Stickleback, row.names = 1,
info.cols = c(2:13), subset.rows = c("Site", "B"))
# Index calculation with jackknife variance estimate
# This can take time for big data sets
Eresult = Eindex(GutContents_SiteB, index = "saramaki", jackknife = TRUE)
rm(list=ls(all=TRUE))
}

\keyword{univar}
\keyword{models}
