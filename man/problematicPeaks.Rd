\name{problematicPeaks}
\alias{problematicPeaks}
\title{Identify intense peaks (in a list of unmatched peaks)}
\usage{
  problematicPeaks(peaks_unmatched, peaks_matched, mode =
    "pH")
}
\arguments{
  \item{peaks_unmatched}{Table of unmatched peaks, with at
  least \code{cpdID, scan, mzFound, int}.}

  \item{peaks_matched}{Table of matched peaks (used for
  base peak reference), with at least \code{cpdID, scan,
  int}.}

  \item{mode}{Processing mode (\code{"pH", "pNa"} etc.)}
}
\value{
  A filtered table with the potentially problematic peaks,
  including the precursor mass and MSMS base peak intensity
  (\code{aMax}) for reference.
}
\description{
  Finds a list of peaks in spectra with a high relative
  intensity (>10% and 1e4, or >1% and 1e5) to write a list
  of peaks which must be manually checked.  Peaks orbiting
  around the parent peak mass (calculated from the compound
  ID), which are very likely co-isolated substances, are
  ignored.
}
\examples{
\dontrun{
# As used in the workflow:
    fp_rean <-  problematicPeaks(
                peaksNoformula,
                specs$peaksMatched,
                mode)
}
}
\author{
  Michael Stravs
}
\seealso{
  \code{\link{msmsWorkflow}}
}

