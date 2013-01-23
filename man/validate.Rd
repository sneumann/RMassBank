\name{validate}
\alias{validate}
\title{Validates a record}
\usage{
  validate(File)
}
\arguments{
  \item{File}{The filepath to the plaintext-records}
}
\value{
  None
}
\description{
  Validates a plaintext-MassBankrecord 
}
\details{
  Reading Annotations does not work at the moment
}
\examples{
validate(c("filepath_to_records/RC00001.txt","filepath_to_records/RC00002.txt"))

}
\author{
  Erik Mueller
}
\seealso{
  \code{\link{parseMassBank}}
}

