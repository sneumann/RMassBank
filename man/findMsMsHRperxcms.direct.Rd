\name{findMsMsHRperxcms.direct}
\alias{findMsMsHRperxcms.direct}
\title{Read in mz-files using XCMS}
\usage{
  findMsMsHRperxcms.direct(fileName, cpdID, mode,
    findPeaksArgs)
}
\arguments{
  \item{fileName}{The path to the mz-file that should be
  read}

  \item{cpdID}{The compoundID of the compound that has been
  used for the file}

  \item{mode}{The ionization mode that has been used for
  the spectrum represented by the peaklist}

  \item{findPeaksArgs}{A list of arguments that will be
  handed to the xcms-method findPeaks via do.call}
}
\value{
  The \code{msmsWorkspace} with the additional peaklist
  added to the right spectrum
}
\description{
  Picks peaks from mz-files and returns the pseudospectra
  that CAMERA creates with the help of XCMS
}
\examples{
\dontrun{
		fileList <- list.files(system.file("XCMSinput", package = "RMassBank"), "Glucolesquerellin", full.names=TRUE)[3]
		loadList(system.file("XCMSinput/compoundList.csv",package="RMassBank"))
     psp <- findMsMsHRperxcms.direct(fileList,2184)
}
}
\author{
  Erik Mueller
}
\seealso{
  \code{\link{msmsWorkflow}}
}

