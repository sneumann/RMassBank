\name{getPcId}
\alias{getPcId}
\title{Search Pubchem CID}
\usage{
  getPcId(search)
}
\arguments{
  \item{search}{The search term.}
}
\value{
  The PubChem CID (in string type).
}
\description{
  Retrieves PubChem CIDs for a search term.
}
\details{
  Only the first result is returned currently. \bold{The
  function should be regarded as experimental and has not
  thoroughly been tested.}
}
\examples{
# Benzene (again):
getPcId("benzene")
}
\author{
  Michael Stravs
}
\references{
  PubChem search: \url{http://pubchem.ncbi.nlm.nih.gov/}

  Entrez E-utilities:
  \url{http://www.ncbi.nlm.nih.gov/books/NBK25500/}
}
\seealso{
  \code{\link{getCtsRecord}}, \code{\link{getCactus}}
}

