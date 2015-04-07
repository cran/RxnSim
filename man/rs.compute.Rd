\name{rs.compute}
\alias{rs.compute}
\alias{rs.compute.list}
\alias{rs.compute.sim.matrix}
\alias{rs.compute.DB}

\title{Computes Similarity of Reactions}

\description{
Computes similarity between two (or more) input reactions.\cr
\code{rs.compute} computes similarity of two reactions.\cr
\code{rs.compute.list} computes similarity of two lists of reactions.\cr
\code{rs.compute.sim.matrix} computes similarity of reactions in a list.\cr
\code{rs.compute.DB} computes similarity of a reaction against a database (parsed from text file).\cr
}

\usage{
rs.compute (rxnA, rxnB, format = 'rsmi', standardize = T, explicitH = F, reversible = T, 
            algo = 'msim', sim.method = 'tanimoto', fp.type = 'extended', fp.mode = 'bit',
            fp.depth = 6, fp.size = 1024, verbose = F, fpCached = F)
rs.compute.list (rxnA, rxnB, format = 'rsmi', standardize = T, explicitH = F, reversible
            = T, algo = 'msim', sim.method = 'tanimoto', fp.type = 'extended', fp.mode 
            = 'bit',fp.depth = 6, fp.size = 1024, clearCache = T)
rs.compute.sim.matrix (rxnA, format = 'rsmi', standardize = T, explicitH = F, reversible
            = T, algo = 'msim', sim.method = 'tanimoto', fp.type = 'extended', fp.mode
            = 'bit', fp.depth = 6, fp.size = 1024, clearCache = T)
rs.compute.DB (rxnA, DB, format = 'rsmi', ecrange = '*', reversible = T, algo = 'msim',
            sim.method = 'tanimoto', sort = T, fpCached = F)
}

\arguments{
  \item{rxnA}{input reaction in RSMI format or name (with path) of MDL RXN file. \code{rs.compute.list} and \code{rs.compute.sim.matrix} accept list of reactions as input.}
  \item{rxnB}{input reaction in RSMI format or name (with path) of MDL RXN file. \code{rs.compute.list} accepts list of reactions as input.}
  \item{DB}{parsed database object as returned by \code{\link[RxnSim:rs.makeDB]{rs.makeDB}}.}
  \item{format}{specifies format of input reaction(s). Reaction(s) can be provided in one of following formats: 'RSMI' (default) or 'RXN'.}
  \item{ecrange}{EC number(s) search pattern while comparing against reaction DB. * is used as wildcard. E.g., 1.2.1.* will restricted search to all reactions with EC numbers starting with 1.2.1.- .}
  \item{standardize}{suppresses all explicit hydrogen if set as \code{TRUE} (default).}
  \item{explicitH}{converts all implicit hydrogen to explicit if set as \code{TRUE}. It is set as \code{FALSE} by default.}
  \item{reversible}{boolean that indicates reversibility of input reaction(s). If set as \code{TRUE} (default), reaction(s) are aligned by comparing them in forward direction and by reversing one of them to compute maximum similarity value.}
  \item{algo}{reaction similarity algorithm to be used. One of following algorithms can be used: \code{'msim'} (default), \code{'msim_max'}, \code{'rsim'} and \code{'rsim2'}. See description for the details of the algorithms.}
  \item{sim.method}{similarity metric to be used to evaluate reaction similarity. Allowed types include:\cr \code{'simple'}, \code{'jaccard'}, \code{'tanimoto'} (default), \code{'russelrao'}, \code{'dice'}, \code{'rodgerstanimoto'}, \code{'achiai'}, \code{'cosine'}, \code{'kulczynski2'}, \code{'mt'}, \code{'baroniurbanibuser'}, \code{'tversky'}, \code{'robust'}, \code{'hamann'}, \code{'pearson'}, \code{'yule'}, \code{'mcconnaughey'}, \code{'simpson'}, \code{'jaccard-count'} and \code{'tanimoto-count'}.}
  \item{fp.type}{fingerprint type to use. Allowed types include:\cr \code{'standard'}, \code{'extended'} (default), \code{'graph'}, \code{'estate'}, \code{'hybridization'}, \code{'maccs'}, \code{'pubchem'}, \code{'kr'}, \code{'shortestpath'}, \code{'signature'} and \code{'circular'}.}
  \item{fp.mode}{fingerprint mode to be used. It can either be set to \code{'bit'} (default) or \code{'count'}.}
  \item{fp.depth}{search depth for fingerprint construction. This argument is ignored for \code{'pubchem'}, \code{'maccs'}, \code{'kr'} and \code{'estate'} fingerprints.}
  \item{fp.size}{length of the fingerprint bit string. This argument is ignored for the \code{'pubchem'}, \code{'maccs'}, \code{'kr'}, \code{'estate'}, \code{'circular'} (count mode) and \code{'signature'} fingerprints.}
  \item{verbose}{boolean that enables display of detailed molecule pairing and reaction alignment (and respective similarity values). The argument is ignored for \code{'rsim2'} algorithm.}
  \item{sort}{boolean than enables \code{rs.compute.DB} to return data frame sorted based upon decreasing value of similarities.}
  \item{fpCached}{boolean that enables fingerprint caching. It is set to \code{FALSE} by default.}
  \item{clearCache}{boolean that resets the cache before (and after) processing reaction lists. It is set to \code{TRUE} by default. Cache can also be explicitly cleared using \code{rs.clearCache}.}
}

\value{
  \item{\code{rs.compute}}{returns a similarity value.}
  \item{\code{rs.compute.list}}{returns a \eqn{m \times n} matrix of similarity values. \eqn{m} and \eqn{n} are the length of two input lists respectively.}
  \item{\code{rs.compute.sim.matrix}}{returns a \eqn{m \times m} symmetric matrix of similarity values. \eqn{m} is the length of the input list.}
  \item{\code{rs.compute.DB}}{returns a data frame.}
}

\details{
\pkg{RxnSim} implements four algorithms to compute reaction similarity, namely \code{msim}, \code{msim_max}, \code{rsim} and \code{rsim2}.
\describe{
  \item{\code{msim}}{is based on individual similarities of molecules in two reactions. First, each reactant (product) of a reaction is paired with an equivalent (similar) reactant (product) of the other reaction based on pairwise similarity values using hierarchical grouping. A \code{0} similarity value is assigned to each unpaired molecule. Reaction similarity is then computed by averaging the similarity values for each pair of equivalent molecule(s) and unpaired molecule(s). Molecule equivalences computed can be reviewed using \code{verbose} mode in \code{rs.compute}.}
  \item{\code{msim_max}}{reaction similarity is computed in the same way as described for \code{msim} except that the unpaired molecules are not used for computing average.}
  \item{\code{rsim}}{is based on cumulative features of reactant(s) and product(s) of two reactions. Each reaction is represented by two fingerprints, one each for the reactants and another for products. Reaction similarity is computed by averaging similarity values obtained by comparing reactants fingerprint and products fingerprints.}
  \item{\code{rsim2}}{is based on cumulative features of all molecules in a reaction forming a reaction fingerprint. Reaction similarity is computed based on the reaction fingerprints of two reactions.}
}

For reversible reactions (\code{reversible = TRUE}), apart from comparing reactions in the forward direction they are also compared by reversing one of the reactions. The greater of the two similarity values is reported.

\strong{Fingerprint Caching}\cr
\code{rs.compute} and \code{rs.compute.DB} functions can use fingerprint caching. If \code{fpCached} is set as \code{TRUE}, cache is queried first before generating fingerprints. Any new fingerprint generated is stored in the cache. Setting \code{fpCached = FALSE} makes no change to cache. Cache can be cleared by calling \code{\link[RxnSim:rs.clearCache]{rs.clearCache}}.

\code{rs.compute.list} and \code{rs.compute.sim.matrix} functions internally use caching. To ensure consistency of fingerprints, \code{\link[RxnSim:rs.clearCache]{rs.clearCache}} is called internally. Use \code{clearCache = FALSE} to override this behaviour; it will use current state of cache and add new fingerprints to it.

Same cache is used for all functions.

\strong{Similarity metric included in \pkg{RxnSim}.} These metric (except \code{jaccard-count} and \code{tanimoto-count}) are derived from \pkg{\link[fingerprint:distance]{fingerprint pacakge}}.
\tabular{lll}{
  \bold{ID} \tab \bold{Name} \tab \bold{Remarks}\cr
  \code{simple} \tab Sokal & Michener \tab bit \cr
  \code{jaccard} \tab Jaccard \tab bit \cr
  \code{tanimoto} \tab Tanimoto (bit) \tab bit and count \cr
  \code{jaccard-count} \tab Jaccard (count) \tab count \cr
  \code{tanimoto-count} \tab Tanimoto (count) \tab count ^\cr
  \code{dice} \tab Dice (bit)\tab bit and count \cr
  \code{russelrao} \tab Russel And Rao \tab bit \cr
  \code{rodgerstanimoto} \tab Roger And Tanimoto \tab bit \cr
  \code{achiai} \tab Ochiai \tab bit \cr
  \code{cosine} \tab Cosine \tab bit \cr
  \code{kulczynski2} \tab Kulczynski 2 \tab bit \cr
  \code{mt} \tab Modified Tanimoto \tab bit \cr
  \code{baroniurbanibuser} \tab Baroni-Urbani/Buser \tab bit \cr
  \code{robust} \tab Robust (bit)\tab bit and count \cr
  \code{tversky} \tab Tversky* \tab bit \cr
  \code{hamann} \tab Hamann \tab bit \cr
  \code{pearson} \tab Pearson \tab bit \cr
  \code{yule} \tab Yule \tab bit \cr
  \code{mcconnaughey} \tab McConnaughey \tab bit \cr
  \code{simpson} \tab Simpson \tab bit \cr
}
*Tversky coefficients can be specified by combining them into a vector, e.g., \code{c('tversky', a, b)}.

\code{tanimoto} (bit), \code{dice} (bit) and \code{robust} (bit) compute similarity of feature vectors (count mode) by translating them to equivalent fingerprint vectors. Default similarity metric used is \code{tanimoto}.\cr\cr
\strong{List of fingerprints included in \pkg{RxnSim}}. These are derived from \link[rcdk:get.fingerprint]{rCDK package}.
\tabular{lll}{
  \bold{ID} \tab \bold{Name of the Fingerprint} \tab \bold{Mode}\cr
  \code{standard} \tab Standard \tab bit \cr
  \code{extended} \tab Extended \tab bit \cr
  \code{estate} \tab EState \tab bit \cr
  \code{graph} \tab Graphonly \tab bit \cr
  \code{hybridization} \tab Hybridization \tab bit \cr
  \code{maccs} \tab MACCS \tab bit \cr
  \code{pubchem} \tab Pubchem \tab bit \cr
  \code{kr} \tab Klekota-Roth \tab bit \cr
  \code{shortestpath} \tab Shortestpath \tab bit \cr
  \code{signature} \tab Signature \tab count \cr
  \code{circular} \tab Circular \tab bit and count \cr
}
}

\examples{
\dontshow{
rct1 <- 'C(=O)C1(=CC=CC=C1).O=O.[OH2]>>[H+].C(C1(C=CC=CC=1))([O-])=O.OO'
rct2 <- 'C(C=CC1(=CC=CC=C1))=O.O=O.[OH2]>>[H+].C(=O)([O-])C=CC1(=CC=CC=C1).OO'
}
\dontrun{Reaction similarity using \'msim\' algorithm}
rs.compute(rct1, rct2, verbose = TRUE)
}

\author{
Varun Giri \email{varungiri@gmail.com}
}

\references{
^ Carbonell, P., Planson, A-G., Fichera, D., & Faulon J-L. (2011) A retrosynthetic biology approach to metabolic pathway design for therapeutic production. \emph{BMC Systems Biology}, \bold{5}:122.
}

\seealso{
\code{\link[RxnSim:rs.makeDB]{rs.makeDB}}, \code{\link[RxnSim:rs.clearCache]{rs.clearCache}}, \code{\link[RxnSim:ms.compute]{ms.compute}}
}

\note{
While using fingerprint caching (by setting \code{fpCached = T} in \code{rs.compute} and \code{rs.compute.DB} or \code{clearCache = F} in \code{rs.compute.list} and \code{rs.compute.sim.matrix}), ensure that the fingerprints are generated using same parameters values (\code{fp.type}, \code{fp.mode}, \code{fp.depth} and \code{fp.size}). To reset cache, call \code{\link[RxnSim:rs.clearCache]{rs.clearCache}}.

\code{rs.compute.DB} uses same parameter values for creating fingerprint as used for (and stored with) DB object (created using \code{\link[RxnSim:rs.makeDB]{rs.makeDB}}) passed as argument.
}