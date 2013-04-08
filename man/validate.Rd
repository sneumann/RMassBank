\name{validate}
\alias{validate}
\title{Validate a record against a set of RUnit tests}
\usage{
  validate(path)
}
\arguments{
  \item{path}{The filepath to a single record, or a
  directory to search recursively}
}
\description{
  Validates a plaintext-MassBankrecord, or recursively all
  records below a directory. The Unit Tests to be used are
  installed in RMassBank/
}
\examples{
\dontrun{
validate("/tmp/MassBank/OpenData/record/")
}
}

