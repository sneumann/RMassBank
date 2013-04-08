\name{addMB}
\alias{addMB}
\title{}
\usage{
  addMB(w, cpdID, fileName, mode)
}
\arguments{
  \item{w}{The \code{msmsWorkspace} to load the peaks into.}

  \item{cpdID}{The compound ID of the substance that has been used for the record}
  
  \item{fileName}{Filename of the MB-record that should be read}
  
  \item{mode}{The ionization mode which was used} 
}
\value{
  The \code{msmsWorkspace} with the added peaks from the record.
}
\description{
  Loads the peaks from a MassBank-record into the supplied msmsWorkspace.
}

\examples{
\dontrun{addMB("PB000001.txt")}
}
\author{
  Erik Mueller
}
\seealso{
  \code{\link{msmsWorkflow}}
}