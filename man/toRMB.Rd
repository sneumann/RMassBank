\name{toRMB}
\alias{toRMB}
\title{Converts XCMS pseudospectra to RMassBank spectrum format}
\usage{
  toRMB(msmsXCMSspecs, cpdID, mode, MS1spec)
}
\arguments{
  \item{msmsXCMSspecs}{A list of pseudospectra from XCMS.}
  
  \item{cpdID}{The compound ID of the substance that has been used to make the pseudospectrum}
  
  \item{mode}{The ionization mode which was used with the substance.}
  
  \item{MS1spec}{The MS1-spectrum in pseudospectrum format}
  
}
\value{
  The converted spectrum in the format that is readable for RMassBank-functions.
}
\description{
  This function converts the format that XCMS pseudospectra have to the format that RMassBank uses for its spectra.
}
\details{
  
}
\examples{
\dontrun{
  toRMB(pseudospectrum)
}
}
\author{
  Erik Mueller
}
\seealso{
  \code{\link{msmsWorkflow}}
}

