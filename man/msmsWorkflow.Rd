\name{msmsWorkflow}
\alias{msmsWorkflow}
\title{RMassBank mass spectrometry pipeline}
\usage{
  msmsWorkflow(w, mode="pH", steps=c(1:8), confirmMode =
    FALSE, newRecalibration = TRUE, useRtLimit = TRUE,
    archivename=NA, readMethod = "mzR", precursorscan.cf =
    FALSE, settings = getOption("RMassBank"), analyzeMethod
    = "formula", progressbar = "progressBarHook")
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

  \item{readMethod}{Several methods are available to get
  peak lists from the input files.  Currently supported are
  "mzR" and "peaklist".  "mzR" reads MS/MS raw data. An
  alternative raw data reader "xcms" is not yet released.
  "peaklist" reads a CSV with two columns and the column
  header "mz", "int".}

  \item{precursorscan.cf}{Whether to fill precursor scan
  entries (cf = carry forward). To be used with files which
  for some reasons do not contain precursor scan IDs in the
  mzML, e.g. AB Sciex converted files.}

  \item{settings}{Options to be used for processing.
  Defaults to the options loaded via
  \code{\link{loadRmbSettings}} et al. Refer to there for
  specific settings.}

  \item{analyzeMethod}{The "method" parameter to pass to
  \code{\link{analyzeMsMs}}.}

  \item{progressbar}{The progress bar callback to use. Only
  needed for specialized applications.  Cf. the
  documentation of \code{\link{progressBarHook}} for
  usage.}
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

