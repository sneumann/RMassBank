\name{analyzeMsMs}
\alias{analyzeMsMs}
\alias{analyzeMsMs.formula}
\alias{analyzeMsMs.intensity}
\title{Analyze MSMS spectra}
\usage{
  analyzeMsMs(msmsPeaks, mode="pH", detail=FALSE,
    run="preliminary", filterSettings =
    getOption("RMassBank")$filterSettings, spectraList =
    getOption("RMassBank")$spectraList, method="formula")

  analyzeMsMs.formula(msmsPeaks, mode="pH", detail=FALSE,
    run="preliminary", filterSettings =
    getOption("RMassBank")$filterSettings, spectraList =
    getOption("RMassBank")$spectraList)

  analyzeMsMs.intensity(msmsPeaks, mode="pH", detail=FALSE,
    run="preliminary", filterSettings =
    getOption("RMassBank")$filterSettings, spectraList =
    getOption("RMassBank")$spectraList)
}
\arguments{
  \item{msmsPeaks}{A group of parent spectrum and
  data-dependent MSMS spectra as returned from
  \code{\link{findMsMsHR}} (refer to the corresponding
  documentation for the precise format specifications).}

  \item{mode}{Specifies the processing mode, i.e. which
  molecule species the spectra contain. \code{\var{pH}}
  (positive H) specifies [M+H]+, \code{\var{pNa}} specifies
  [M+Na]+, \code{\var{pM}} specifies [M]+, \code{\var{mH}}
  and \code{\var{mNa}} specify [M-H]- and [M-Na]-,
  respectively. (I apologize for the naming of
  \code{\var{pH}} which has absolutely nothing to do with
  chemical \emph{pH} values.)}

  \item{detail}{Whether detailed return information should
  be provided (defaults to \code{FALSE}). See below.}

  \item{run}{\code{"preliminary"} or \code{"recalibrated"}.
  In the \code{preliminary} run, mass tolerance is set to
  10 ppm (above m/z 120) and 15 ppm (below m/z 120), the
  default intensity cutoff is $10^4$ for positive mode (no
  default cutoff in negative mode), and the column
  \code{"mz"} from the spectra is used as data source.  In
  the \code{recalibrated} run, the mass tolerance is set to
  5 ppm over the whole mass range, the default cutoff is 0
  and the column \code{"mzRecal"} is used as source for the
  m/z values. Defaults to \code{"preliminary"}.}

  \item{filterSettings}{Settings for the filter parameters,
  by default loaded from the RMassBank settings set with
  e.g. \code{\link{loadRmbSettings}}. Must contain:
  \itemize{ \item \code{ppmHighMass}, allowed ppm deviation
  before recalibration for high mass range \item
  \code{ppmLowMass}, allowed ppm deviation before
  recalibration for low mass range \item
  \code{massRangeDivision}, division point between high and
  low mass range (before recalibration) \item
  \code{ppmFine}, allowed ppm deviation overall after
  recalibration \item \code{prelimCut}, intensity cutoff
  for peaks in preliminary run \item \code{prelimCutRatio},
  relative intensity cutoff for peaks in preliminary run,
  e.g. 0.01 = 1% \item \code{fineCut}, intensity cutoff for
  peaks in second run \item \code{fineCutRatio}, relative
  intensity cutoff for peaks in second run \item
  \code{specOkLimit}, minimum intensity of base peak for
  spectrum to be accepted for processing \item
  \code{dbeMinLimit}, minimum double bond equivalent for
  accepted molecular subformula.  \item
  \code{satelliteMzLimit}, for satellite peak filtering
  (\code{\link{filterPeakSatellites}}: mass window to use
  for satellite removal \item \code{satelliteIntLimit}, the
  relative intensity below which to discard "satellites".
  (refer to \code{\link{filterPeakSatellites}}).  }}

  \item{spectraList}{The list of MS/MS spectra present in
  each data block. As also defined in the settings file.}

  \item{method}{Selects which function to actually use for
  data evaluation. The default "formula" runs a full
  analysis via formula assignment to fragment peaks. The
  alternative setting "intensity" calls a "mock"
  implementation which circumvents formula assignment and
  filters peaks purely based on intensity cutoffs and the
  satellite filtering. (In this case, the ppm and dbe
  related settings in filterSettings are ignored.)}
}
\value{
  \item{list("foundOK")}{ Boolean. Whether or not child
  spectra are present for this compound (inherited from
  \code{msmsdata}).} \item{list("mzrange")}{ The maximum
  m/z range over all child spectra.} \item{list("id")}{ The
  compound ID (inherited from \code{msmsdata})}
  \item{list("mode")}{ processing mode}
  \item{list("parentHeader")}{ Parent spectrum header data
  (ex \code{msmsdata})} \item{list("parentMs")}{ Parent
  spectrum (ex \code{msmsdata}) in matrix format}
  \item{list("msmsdata")}{ Analysis results for all child
  spectra: \itemize{ \item\code{specOK} Boolean. Whether or
  not the spectrum contains any useful peaks. If
  \code{specOK = FALSE}, all other information (except scan
  info and compound ID) may be missing!  \item\code{parent}
  Parent mass and formula in a one-row data frame format.
  Currently rather obsolete, originally contained data from
  MolgenMsMs results.  \item \code{childFilt} Annotated
  peaks of the MSMS spectrum (after filtering by accuracy)
  \item \code{childRaw} Raw (\code{mz, int}) spectrum
  before any treatment. (With recalibrated data, this is
  (\code{mz, int, mzRecal}).  } For \code{detail = TRUE},
  additionally: \itemize{ \item\code{childRawLow} Peaks cut
  away because of low (absolute or relative) intensity
  \item\code{childRawSatellite} Peaks cut away
  as"satellites" \item\code{childRawOK} Peaks after cutting
  away low/satellite peaks. Used for further analysis steps
  \item\code{child} Annotated peaks of the MSMS spectrum
  before filtering by accuracy \item \code{childBad}
  Annotated peaks of the MSMS spectrum which didn't pass
  the accuracy threshold \item\code{childUnmatched} Peaks
  of the MSMS spectrum with no annotated formula }}
}
\description{
  Analyzes MSMS spectra of a compound by fitting formulas
  to each subpeak.
}
\details{
  The analysis function uses Rcdk. Note that in this step,
  \emph{satellite peaks} are removed by a simple heuristic
  rule (refer to the documentation of
  \code{\link{filterPeakSatellites}} for details.)
}
\examples{
\dontrun{analyzed <- analyzeMsMs(spec, "pH", TRUE)}
}
\author{
  Michael Stravs
}
\seealso{
  \code{\link{msmsWorkflow}},
  \code{\link{filterLowaccResults}},
  \code{\link{filterPeakSatellites}},
  \code{\link{reanalyzeFailpeaks}}
}

