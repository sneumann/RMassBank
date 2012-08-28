\name{getCtsRecord}
\alias{getCtsRecord}
\title{Retrieve information from CTS}
\usage{
  getCtsRecord(key, from = "inchikey", to = c("cas",
    "hmdb", "kegg", "sid", "chebi", "inchi", "lipidmap",
    "smiles", "cid", "inchikey", "mass", "formula",
    "iupac", "names"))
}
\arguments{
  \item{key}{The search term (or key).}

  \item{from}{The format of the key. Allowed are
  \code{"cas", "hmdb", "kegg", "sid", "chebi", "inchi",
  "lipidmap", "smiles", "cid", "inchikey", "mass",
  "formula", "iupac", "name"}.}

  \item{to}{The list of result types which should be
  returned. Allowed are \code{"cas", "hmdb", "kegg", "sid",
  "chebi", "inchi", "lipidmap", "smiles", "cid",
  "inchikey", "mass", "formula", "iupac", "name"}.}
}
\value{
  Returns a named list with the values of the results. The
  list item \code{"names"} is a matrix with columns
  \code{"name", "score"}, with \code{score} being an
  indicator of the reliability of the name assignment.
}
\description{
  Retrieves chemical information about a compound from
  Chemical Translation Service (CTS) from a known
  identifier.
}
\note{
  The return values are not 100% reliable, e.g. a known bug
  returns "ChEBI" for the \code{chebi} entry instead of the
  actual ChEBI code in some instances.
}
\examples{
getCtsRecord("benzene", "name")
}
\author{
  Michael Stravs
}
\references{
  Chemical Translation Service:
  \url{http://uranus.fiehnlab.ucdavis.edu:8080/cts/homePage}
}
\seealso{
  \code{\link{getCactus}},\code{\link{getPcId}}
}

