\name{aggregateSpectra}
\alias{aggregateSpectra}
\title{Aggregate analyzed spectra}
\usage{
  aggregateSpectra(spec, addIncomplete = FALSE)
}
\arguments{
  \item{spec}{The set of spectra to aggregate}

  \item{addIncomplete}{Whether or not the peaks from
  incomplete files (files for which less than the maximal
  number of spectra are present)}
}
\value{
  \item{foundOK }{ A numeric vector with the compound IDs
  of all files for which spectra were found.
  \code{names(foundOK)} are the filenames.} \item{foundFail
  }{ A numeric vector with the compound IDs of all files
  for which no spectra were found. \code{names(foundOK)}
  are the filenames.} \item{spectraFound }{ A numeric
  vector indicated the number of found spectra per
  compound} \item{specFound }{ A list of processed spectral
  data for all compounds with at least 1 found spectrum, as
  returned by \code{\link{analyzeMsMs}}.} \item{specEmpty
  }{ A list of (not-really-)processed spectral data for
  compounds without spectra.} \item{specComplete }{ A list
  of processed spectral data for all compounds with the
  full spectrum count (i.e.
  \code{length(getOption("RMassBank")$spectraList)}
  spectra.)  As such, \code{specComplete} is a subset of
  \code{specFound}.} \item{specIncomplete}{ A list of
  processed spectral data for all compounds with incomplete
  spectrum count. The complement to \code{specComplete}.}
  \item{peaksMatched }{ A dataframe of all peaks with a
  matched formula, which survived the elimination criteria.
  } \item{peaksUnmatched }{ A dataframe of all peaks
  without a matched formula, or with a formula which failed
  the filter criteria.}
}
\description{
  Groups an array of analyzed spectra and creates
  aggregated peak tables
}
\details{
  \code{\var{addIncomplete}} is relevant for recalibration.
  For recalibration, we want to use only high-confidence
  peaks, therefore we set \code{\var{addIncomplete}} to
  \code{FALSE}. When we want to generate a peak list for
  actually generating MassBank records, we want to include
  all peaks into the peak tables.
}
\examples{
## As used in the workflow:
\dontrun{\%
	analyzedRcSpecs <- lapply(recalibratedSpecs, function(spec)
		analyzeMsMs(spec, mode="pH", detail=TRUE, run="recalibrated", cut=0, cut_ratio=0 ) )
	aggregatedSpecs <- aggregateSpectra(analyzedSpecs)
}
}
\author{
  Michael Stravs
}
\seealso{
  \code{\link{msmsWorkflow}}, \code{\link{analyzeMsMs}}
}

