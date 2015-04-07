\name{rs.mask}
\alias{rs.mask}
\alias{ms.mask}

\title{Masks a Sub-structure in Input Molecule or Reaction}

\description{
Replaces a sub-structure, provided as SMILES or SMARTS, with the given mask.\cr
\code{rs.mask} masks input sub-structure in reaction.\cr
\code{ms.mask} masks input sub-structure in molecule.\cr
}

\usage{
rs.mask (substructure, mask, reaction, format = 'rsmi', standardize = T,
         explicitH = F, recursive = F)
ms.mask (substructure, mask, molecule, format = 'smiles', standardize = T,
         explicitH = F, recursive = F)
}

\arguments{
  \item{substructure}{SMILES or SMARTS of the structure to be searched and masked.}
  \item{mask}{SMILES of structure to be used as mask.}
  \item{reaction}{Input reaction to be processed.}
  \item{molecule}{Input molecule to be processed.}
  \item{format}{specifies format of input reaction/molecule. It can be one of following for a reaction: 'RSMI' or 'RXN'; for a molecule: 'SMILES' or 'MOL'.}
  \item{standardize}{suppresses all explicit hydrogen if set as \code{TRUE} (default).}
  \item{explicitH}{converts all implicit hydrogen to explicit if set as \code{TRUE}. It is set as \code{FALSE} by default.}
  \item{recursive}{if TRUE, all the occurrences of input sub-structure are replaced.}
}

\value{
  Returns SMILES with mask.
}

\details{
The sub-structure is searched in input reaction/molecule and replaced with the mask. All the bonds between identified sub-structure and the remaining atoms are mapped to the mask. If mask contains more than one atom, all the bonds are connected to the last atom in mask. By default, the first identified sub-structure is replaced. To replace all occurrences, recursive should be set to \code{TRUE}. Valence is not checked for the mask atom and the final structure.
}

\examples{
ms.mask('OP(=O)O', '[Cs]', 'O=P(O)(O)OP(=O)(O)OP(=O)(O)OCC3OC(n2cnc1c(ncnc12)N)C(O)C3O')
}

\author{
Varun Giri \email{varungiri@gmail.com}
}

\seealso{
\code{\link[RxnSim:rs.makeDB]{rs.makeDB}}
}

\note{
Aromatic form of SMILES of the query sub-structure should be used for masking aromatic structures. Automatic aromaticity perception is not done on query structures.
}