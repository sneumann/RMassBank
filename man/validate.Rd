\name{validate}
\alias{validate}
\title{Validate MassBank records with a set of Unit tests}
\usage{
  validate(path)
}
\arguments{
  \item{path}{The filepath to a single record, or a
  directory to search recursively}
}
\description{
  Validates a plain text MassBank record, or recursively
  all records within a directory. The Unit Tests to be used
  are installed in RMassBank/inst/unitTests and currently
  include checks for NAs, peaks versus precursor, precursor
  mz, precursor type, SMILES vs exact mass, total
  intensities and title versus type. The validation report
  is saved as "report.html" in the working directory.
}
\examples{
\dontrun{
validate("/tmp/MassBank/OpenData/record/")
}
}

