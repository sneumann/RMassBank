\name{reanalyzeFailpeaks}
\alias{reanalyzeFailpeak}
\alias{reanalyzeFailpeaks}
\title{Reanalyze unmatched peaks}
\usage{
  reanalyzeFailpeaks(specs, custom_additions, mode,
    filterSettings = getOption("RMassBank")$filterSettings)
    reanalyzeFailpeak(custom_additions, mass, cpdID,
    counter, pb = NULL, mode, filterSettings =
    getOption("RMassBank")$filterSettings)
}
\arguments{
  \item{specs}{An \code{aggregatedRcSpecs} object (after
  the electronic noise was cleared from the unmatched
  peaks).}

  \item{custom_additions}{The allowed additions, e.g.
  "N2O".}

  \item{mode}{Processing mode (\code{"pH", "pNa", "mH"}
  etc.)}

  \item{mass}{(Usually recalibrated) m/z value of the
  peak.}

  \item{cpdID}{Compound ID of this spectrum.}

  \item{counter}{Current peak index (used exclusively for
  the progress indicator)}

  \item{pb}{A txtProgressBar object to display progress on.
  No progress is displayed if NULL.}

  \item{filterSettings}{Settings for filtering data. Refer
  to\code{\link{analyzeMsMs}} for settings.}
}
\value{
  The returning list contains two tables:
  \item{peaksReanalyzed}{All reanalyzed peaks with or
  without matching formula.}
  \item{peaksMatchedReanalysis}{Only the peaks with a
  matched reanalysis formula.}

  It would be good to merge the analysis functions of
  \code{analyzeMsMs} with the one used here, to simplify
  code changes.
}
\description{
  Reanalysis of peaks with no matching molecular formula by
  allowing additional elements (e.g. "N2O").
}
\details{
  \code{reanalyzeFailpeaks} examines the
  \code{unmatchedPeaksC} table in \code{specs} and sends
  every peak through \code{reanalyzeFailpeak}.
}
\examples{
## As used in the workflow:
\dontrun{
	reanalyzedRcSpecs <- reanalyzeFailpeaks(aggregatedRcSpecs, custom_additions="N2O", mode="pH")
# A single peak:
reanalyzeFailpeak("N2O", 105.0447, 1234, 1, 1, "pH")
}
}
\author{
  Michael Stravs
}
\seealso{
  \code{\link{analyzeMsMs}}, \code{\link{msmsWorkflow}}
}

