\name{msmsWorkflow}
\alias{msmsWorkflow}
\title{RMassBank mass spectrometry pipeline}
\usage{
  msmsWorkflow(w, mode = "pH", steps = c(1:8), confirmMode
    = FALSE, newRecalibration = TRUE, useRtLimit = TRUE,
    archivename = NA)
}
\arguments{
  \item{w}{A \code{msmsWorkspace} to work with.}

  \item{mode}{\code{"pH", "pNa", "pM", "mH", "mM", "mFA"}
  for different ions ([M+H]+, [M+Na]+, [M]+, [M-H]-, [M]-,
  [M+FA]-).}

  \item{steps}{Which steps of the workflow to process. See
  the vignette \code{vignette("RMassBank")} for details.}

  \item{confirmMode}{Defaults to false (use most intense
  precursor). Value 1 uses the 2nd-most intense precursor
  for a chosen ion (and its data-dependent scans) , etc.}

  \item{newRecalibration}{Whether to generate a new
  recalibration curve (\code{TRUE}, default) or to reuse
  the currently stored curve (\code{FALSE}, useful e.g. for
  adduct-processing runs.)}

  \item{useRtLimit}{Whether to enforce the given retention
  time window.}

  \item{archivename}{The prefix under which to store the
  analyzed result files.}
}
\value{
  The processed \code{msmsWorkspace}.
}
\description{
  Extracts and processes spectra from a specified file
  list, according to loaded options and given parameters.
}
\details{
  The filenames of the raw LC-MS runs are read from the
  array \code{files} in the global enviroment. See the
  vignette \code{vignette("RMassBank")} for further details
  about the workflow.
}
\author{
  Michael Stravs, Eawag <michael.stravs@eawag.ch>
}
\seealso{
  \code{\link{msmsWorkspace-class}}
}

