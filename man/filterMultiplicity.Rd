\name{filterMultiplicity}
\alias{filterMultiplicity}
\title{filterMultiplicity}
\usage{
  filterMultiplicity(specs, archivename = NA, mode = "pH",
    recalcBest = TRUE)
}
\arguments{
  \item{specs}{aggregatedSpecs object whose peaks should be
  filtered}

  \item{archivename}{The archive name, used for generation
  of archivename_failpeaks.csv}

  \item{mode}{Mode of ion analysis}

  \item{recalcBest}{Boolean, whether to recalculate the
  formula multiplicity after the first multiplicity
  filtering step. Sometimes, setting this to FALSE can be a
  solution if you have many compounds with e.g. fluorine
  atoms, which often have multiple assigned formulas per
  peak and might occasionally lose peaks because of that.}
}
\value{
  A list object with values: \item{peaksOK}{ Peaks with
  >1-fold formula multiplicity from the "normal" peak
  analysis.  } \item{peaksReanOK}{ Peaks with >1-fold
  formula multiplicity from peak reanalysis.  }
  \item{peaksFiltered}{ All peaks with annotated formula
  multiplicity from first analysis.  }
  \item{peaksFilteredReanalysis}{ All peaks with annotated
  formula multiplicity from peak reanalysis.  }
  \item{peaksProblematic}{ Peaks with high intensity which
  do not match inclusion criteria -> possible false
  negatives. The list will be exported into
  archivename_failpeaks.csv. }
}
\description{
  Multiplicity filtering: Removes peaks which occur only
  once in a n-spectra set.
}
\details{
  This function executes multiplicity filtering for a set
  of spectra using the workhorse function
  \code{\link{filterPeaksMultiplicity}} (see details there)
  and retrieves problematic filtered peaks (peaks which are
  of high intensity but were discarded, because either no
  formula was assigned or it was not present at least 2x),
  using the workhorse function
  \code{\link{problematicPeaks}}. The results are returned
  in a format ready for further processing with
  \code{\link{mbWorkflow}}.
}
\examples{
\dontrun{
    refilteredRcSpecs <- filterMultiplicity(
			reanalyzedRcSpecs, "myarchive", "pH")
}
}
\author{
  Michael Stravs
}
\seealso{
  \code{\link{filterPeaksMultiplicity}},\code{\link{problematicPeaks}}
}

