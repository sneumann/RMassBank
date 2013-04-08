\name{addPeaksManually}
\alias{addPeaksManually}
\title{}
\usage{
  addPeaksManually(w, cpdID, handSpec, mode)
}
\arguments{
  \item{w}{The \code{msmsWorkspace} to load the peaks into.}

  \item{cpdID}{The compound ID of the substance that has been used for the spectra}
  
  \item{handSpec}{Filename of the csv file}
  
  \item{mode}{The ionization mode which was used} 
}
\value{
  The \code{msmsWorkspace} with the added peaks.
}
\description{
  Loads a table with peaks to add to the
  msms-spectra. Required columns are \code{mz, int}.
}

\examples{
\dontrun{addPeaksManually("peaks.csv")}
}
\author{
  Erik Mueller
}
\seealso{
  \code{\link{msmsWorkflow}}
}