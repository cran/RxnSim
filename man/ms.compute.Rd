\name{ms.compute}
\alias{ms.compute}
\alias{ms.compute.sim.matrix}

\title{Computes Similarity of Molecules}

\description{
Computes chemical similarity between two (or more) input molecules.\cr
}

\usage{
ms.compute (molA, molB, format = 'smiles', standardize = T, explicitH = F, sim.method =
            'tanimoto', fp.type = 'extended', fp.mode = 'bit', fp.depth = 6, fp.size = 
            1024, fpCached = F)
ms.compute.sim.matrix (molA, format = 'smiles', standardize = T, explicitH = F,
            sim.method = 'tanimoto', fp.type = 'extended', fp.mode = 'bit', fp.depth = 6,
            fp.size = 1024, clearCache = T)
}

\arguments{
  \item{molA}{input molecule in SMILES format or name (with path) of MDL MOL file. \code{ms.compute.sim.matrix} accepts list of molecules as input.}
  \item{molB}{input molecule in SMILES format or name (with path) of MDL MOL file.}
  \item{format}{specifies format of input molecule(s). Molecule(s) can be provided in one of following formats: \code{'SMILES'} (default) or \code{'MOL'}.}
  \item{standardize}{suppresses all explicit hydrogen if set as \code{TRUE} (default).}
  \item{explicitH}{converts all implicit hydrogen to explicit if set as \code{TRUE}. It is set as \code{FALSE} by default.}
  \item{sim.method}{similarity metric to be used to evaluate molecule similarity. Allowed types include:\cr \code{'simple'}, \code{'jaccard'}, \code{'tanimoto'} (default), \code{'russelrao'}, \code{'dice'}, \code{'rodgerstanimoto'}, \code{'achiai'}, \code{'cosine'}, \code{'kulczynski2'}, \code{'mt'}, \code{'baroniurbanibuser'}, \code{'tversky'}, \code{'robust'}, \code{'hamann'}, \code{'pearson'}, \code{'yule'}, \code{'mcconnaughey'}, \code{'simpson'}, \code{'jaccard-count'} and \code{'tanimoto-count'}.}
  \item{fp.type}{fingerprint type to use. Allowed types include:\cr \code{'standard'}, \code{'extended'} (default), \code{'graph'}, \code{'estate'}, \code{'hybridization'}, \code{'maccs'}, \code{'pubchem'}, \code{'kr'}, \code{'shortestpath'}, \code{'signature'} and \code{'circular'}.}
  \item{fp.mode}{fingerprint mode to be used. It can either be set to \code{'bit'} (default) or \code{'count'}.}
  \item{fp.depth}{search depth for fingerprint construction. This argument is ignored for \code{'pubchem'}, \code{'maccs'}, \code{'kr'} and \code{'estate'} fingerprints.}
  \item{fp.size}{length of the fingerprint bit string. This argument is ignored for \code{'pubchem'}, \code{'maccs'}, \code{'kr'}, \code{'estate'}, \code{'circular'} (count mode) and \code{'signature'} fingerprints.}
  \item{fpCached}{boolean that enables fingerprint caching. It is set to \code{FALSE} by default.}
  \item{clearCache}{boolean that resets the cache before (and after) processing molecule lists.  It is set to \code{TRUE} by default. Cache can also be explicitly cleared by using \code{\link[RxnSim:rs.clearCache]{rs.clearCache}}.}
}

\value{
  Returns similarity value(s).
  \item{\code{ms.compute}}{returns a similarity value.}
  \item{\code{ms.compute.sim.matrix}}{returns a \eqn{m \times m} symmetric matrix of similarity values. \eqn{m} is the length of the input list.}
}

\details{
See \code{\link[RxnSim:rs.compute]{rs.compute}} functions, for details for fingerprints and similarity matrices.
\code{ms.compute} can use fingerprint caching by enabling \code{fpCached} option. \code{ms.compute} and \code{ms.compute.sim.matrix} use same cache as \code{rs.compute} and other functions in the package.
}

\examples{
ms.compute('N', '[H]N([H])[H]', standardize = FALSE)
}

\author{
Varun Giri \email{varungiri@gmail.com}
}

\seealso{
\code{\link[RxnSim:rs.compute]{rs.compute}}, \code{\link[RxnSim:rs.clearCache]{rs.clearCache}}
}

\note{
Fingerprint cache stores fingerprints generated for a molecule index based on its SMILES. When caching is enabled, the fingerprint for a molecule, if present, is retrieved from the cache. The parameters pertaining to fingerprint generation are thus ignored. If the fingerprint for the molecule is not already cached, fingerprint based on the input parameters is generated and stored in the cache.
}