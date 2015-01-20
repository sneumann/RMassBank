\name{getCtsKey}
\alias{getCtsKey}
\title{Convert a single ID to another using CTS.}
\usage{
  getCtsKey(query, from = "Chemical Name", to = "InChIKey")
}
\arguments{
  \item{query}{ID to be converted}

  \item{from}{Type of input ID}

  \item{to}{Desired output ID}
}
\value{
  An unordered array with the resulting converted key(s).
}
\description{
  Convert a single ID to another using CTS.
}
\examples{
\dontrun{
k <- getCtsKey("benzene", "Chemical Name", "InChIKey")
}
}
\author{
  Michele Stravs, Eawag <stravsmi@eawag.ch>
}

