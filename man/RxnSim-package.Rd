\name{RxnSim-package}
\alias{RxnSim-package}
\alias{RxnSim}
\alias{rxnsim}
\docType{package}

\title{Functions to compute chemical reaction and molecular similarity}

\description{
RxnSim provides methods to compute molecular and reaction similarity. It uses rCDK package (R Interface to the CDK Libraries) and fingerprints package for chemoinformatic routines.
}

\details{
RxnSim provides methods to compute chemical similarity between two or more reactions and molecules. Molecular similarity is computed based on structural features. Reaction similarity is a function of similarities of participating molecules. The package provides multiple methods to extract structural features as fingerprints (or feature vectors) and similarity metrics. It additionally provides functionality to mask chemical substructures for weighted similarity computations. It uses rCDK and fingerprint packages for cheminformatics functionality. 

\strong{User functions:}
\describe{
  \item{\code{\link[RxnSim:rs.compute]{rs.compute}}}{computes similarity between two reactions.}
  \item{\code{\link[RxnSim:rs.compute]{rs.compute.list}}}{computes similarity between all pairs of reactions from two lists.}
  \item{\code{\link[RxnSim:rs.compute]{rs.compute.sim.matrix}}}{computes pairwise similarity between all reactions in a list.}
  \item{\code{\link[RxnSim:rs.compute]{rs.compute.DB}}}{computes similarity of a reaction to those in a reaction database (DB) object read from a text file.}
  \item{\code{\link[RxnSim:rs.makeDB]{rs.makeDB}}}{reads a text file containing EC Numbers, Reaction Names and Reaction SMILES and converts it into a reaction DB object.}
  \item{\code{\link[RxnSim:ms.compute]{ms.compute}}}{computes similarity between two molecules.}
  \item{\code{\link[RxnSim:ms.compute.sim.matrix]{ms.compute.sim.matrix}}}{computes pairwise similarity between all molecules in a list.}
  \item{\code{\link[RxnSim:rs.clearCache]{rs.clearCache}}}{clears fingerprint cache.}
  \item{\code{\link[RxnSim:rs.mask]{rs.mask}}}{substitutes given sub-structure in the molecules of a reaction by a user defined mask.}
  \item{\code{\link[RxnSim:ms.mask]{ms.mask}}}{substitutes given sub-structure in a molecule by a user defined mask.}
}
}

\author{
Varun Giri \email{varungiri@gmail.com}

Maintainer: Varun Giri
}

\seealso{
\code{\link[RxnSim:rs.compute]{rs.compute}}, \code{\link[RxnSim:ms.compute]{ms.compute}} 
}

\examples{
\dontrun{Reaction similarity}
rs.compute('CCCO>>CCC=O', 'CC(O)C>>CC(=O)C')

\dontrun{Metabolite similarity}
ms.compute('CCC=O', 'CC(O)C')
}
