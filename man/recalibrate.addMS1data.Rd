\name{recalibrate.addMS1data}
\alias{recalibrate.addMS1data}
\title{Return MS1 peaks to be used for recalibration}
\usage{
  recalibrate.addMS1data(spec,mode="pH",
    recalibrateMS1Window =
    getOption("RMassBank")$recalibrateMS1Window)
}
\arguments{
  \item{spec}{A \code{aggregatedSpecs}-like object.}

  \item{mode}{\code{"pH", "pNa", "pM", "mH", "mM", "mFA"}
  for different ions ([M+H]+, [M+Na]+, [M]+, [M-H]-, [M]-,
  [M+FA]-).}

  \item{recalibrateMS1Window}{Window width to look for MS1
  peaks to recalibrate (in ppm).}
}
\value{
  A dataframe with columns \code{mzFound, formula, mzCalc,
  dppm, dbe, int, dppmBest, formulaCount, good, cpdID,
  scan, parentScan, dppmRc}. However, columns \code{dbe,
  int, formulaCount, good, scan, parentScan} do not contain
  real information and are provided only as fillers.
}
\description{
  Returns the precursor peaks for all MS1 spectra in the
  \code{spec} dataset with annotated formula to be used in
  recalibration.
}
\details{
  For all spectra in \code{spec$specFound}, the precursor
  ion is extracted from the MS1 precursor spectrum. All
  found ions are returned in a data frame with a format
  matching \code{spec$peaksMatched} and therefore suitable
  for \code{rbind}ing to the \code{spec$peaksMatched}
  table. However, only minimal information needed for
  recalibration is returned.
}
\examples{
\dontrun{
# More or less as used in recalibrateSpectra:
		rcdata <- subset(aggregatedSpecs$peaksMatched, formulaCount==1)
		ms1data <- recalibrate.addMS1data(aggregatedSpecs, "pH", 15)
		rcdata <- rbind(rcdata, ms1data)
 # ... continue constructing recalibration curve with rcdata
}
}
\author{
  Michael Stravs, EAWAG <michael.stravs@eawag.ch>
}

