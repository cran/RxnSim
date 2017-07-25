\name{rs.makeDB}
\alias{rs.makeDB}

\title{Converts Text File to Reaction Database}

\description{
Reads and parses input text file containing reaction smiles into reaction database object. The reaction database is used for querying reaction similarity of candidate reactions.
}

\usage{
rs.makeDB (txtFile, header = F, sep = '\t', standardize = T, explicitH = F, fp.type = 
          'extended', fp.mode = 'bit', fp.depth = 6, fp.size = 1024, useMask = F,
          maskStructure, mask, recursive = F)
}

\arguments{
  \item{txtFile}{input file containing EC numbers, reaction name and RSMI. See description for format of input file.}
  \item{header}{boolean to indicate if the input file contains a header. It is set to \code{FALSE} by default.}
  \item{sep}{the field separator character to be used while reading the input file.}
  \item{standardize}{suppresses all explicit hydrogen if set as \code{TRUE} (default).}
  \item{explicitH}{converts all implicit hydrogen to explicit if set as \code{TRUE}. It is set as \code{FALSE} by default.}
  \item{fp.type}{Fingerprint type to use. Allowed types include:\cr \code{'standard'}, \code{'extended'} (default), \code{'graph'}, \code{'estate'}, \code{'hybridization'}, \code{'maccs'}, \code{'pubchem'}, \code{'kr'}, \code{'shortestpath'}, \code{'signature'} and \code{'circular'}.}
  \item{fp.mode}{fingerprint mode to be used. It can either be set to \code{'bit'} (default) or \code{'count'}.}
  \item{fp.depth}{search depth for fingerprint construction. This argument is ignored for \code{'pubchem'}, \code{'maccs'}, \code{'kr'} and \code{'estate'} fingerprints.}
  \item{fp.size}{length of the fingerprint bit string. This argument is ignored for \code{'pubchem'}, \code{'maccs'}, \code{'kr'}, \code{'estate'}, \code{'circular'} (count mode) and \code{'signature'} fingerprints.}
  \item{useMask}{boolean to indicate use of masking. If \code{TRUE}, each reaction is processed to mask given substructure. See \code{\link[RxnSim:rs.mask]{rs.mask}} for details.}
  \item{maskStructure}{SMILES or SMARTS of the structure to be searched and masked.}
  \item{mask}{SMILES of structure to be used as mask.}
  \item{recursive}{if \code{TRUE}, all the occurrences of input substructure are replaced recursively.}
}

\value{
Returns a list, containing parsed input data, reaction fingerprints.
\item{Data}{data frame containing EC Numbers, Reaction Names and RSMI as read from the input file. MaskedRSMI are also included if masking is used.}
\item{FP}{list of molecular fingerprints for each reaction in the input file. These fingerprints are further processed based on the reaction similarity algorithm.}
It also contains the parameter values used for generating fingerprints, viz., \code{standardize}, \code{explicitH}, \code{fp.type},  \code{fp.mode}, \code{fp.depth} and \code{fp.size}.
}

\details{
The parameters used to generate fingerprints are stored in the database object and returned with the parsed data. Same parameter values are used while parsing input reaction in \code{\link[RxnSim:rs.compute.DB]{rs.compute.DB}}.

The input text file should contain following three fields, separated with  \code{TAB} (or any appropriate field separator). A field can be left blank.
\tabular{lll}{
[EC Number] \tab [Reaction Name] \tab [Reaction SMILES (RSMI)]
}

The package comes with a sample reaction database file extracted from Rhea database (Morgat et al., 2015). If no \code{textfile} is provided, default sample database file is used:
\preformatted{
rs.makeDB()
}
A larger dataset containing all reactions from Rhea database (v.83) is also provided with the package.
}

\author{
Varun Giri \email{varungiri@gmail.com}
}

\references{
Morgat, A., Lombardot, T., Axelsen, K., Aimo, L., Niknejad, A., Hyka-Nouspikel, N., Coudert, E., Pozzato, M., Pagni, M., Moretti, S., Rosanoff, S., Onwubiko, J., Bougueleret, L., Xenarios, I., Redaschi, N., Bridge, A. (2017) Updates in Rhea - an expert curated resource of biochemical reactions. \emph{Nucleic Acids Research}, \bold{45}:D415-D418; doi: 10.1093/nar/gkw990
}

\seealso{
\code{\link[RxnSim:rs.compute.DB]{rs.compute.DB}}, \code{\link[RxnSim:rs.mask]{rs.mask}}
}