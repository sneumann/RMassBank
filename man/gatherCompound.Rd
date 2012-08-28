\name{gatherCompound}
\alias{gatherCompound}
\alias{gatherSpectrum}
\title{Compose data block of MassBank record}
\usage{
  gatherCompound(spec, refiltered, additionalPeaks = NULL)

  gatherSpectrum(spec, msmsdata, ac_ms, ac_lc, refiltered,
    additionalPeaks = NULL)
}
\arguments{
  \item{spec}{An object of "analyzedSpectrum" type (i.e.
  contains \code{info}, \code{mzrange}, a list of
  \code{msmsdata}, compound ID, parent MS1, cpd id...)}

  \item{refiltered}{The \code{refilteredRcSpecs} dataset
  which contains our good peaks.  Contains \code{peaksOK},
  \code{peaksReanOK}, \code{peaksFiltered},
  \code{peaksFilteredReanalysis}, \code{peaksProblematic}.
  Currently we use \code{peaksOK} and \code{peaksReanOK} to
  create the spectra.}

  \item{msmsdata}{The \code{msmsdata} sub-object from the
  compound's \code{spec} which is the child scan which is
  currently processed.  Contains \code{childFilt,
  childBad}, scan number, etc. Note that the peaks are
  actually not taken from this list! They were taken from
  \code{msmsdata} initially, but after introduction of the
  refiltration and multiplicity filtering, this was
  changed. Now only the scan information is actually taken
  from \code{msmsdata}.}

  \item{ac_ms,ac_lc}{Information for the
  AC\$MASS_SPECTROMETRY and AC\$CHROMATOGRAPHY fields in
  the MassBank record, created by \code{gatherCompound} and
  then fed into \code{gatherSpectrum}.}

  \item{additionalPeaks}{If present, a table with
  additional peaks to add into the spectra.  As loaded with
  \code{\link{addPeaks}}.}
}
\value{
  \code{gatherCompound} returns a list of tree-like
  MassBank data blocks. \code{gatherSpectrum} returns one
  single MassBank data block or \code{NA} if no useful peak
  is in the spectrum.
}
\description{
  \code{gatherCompound} composes the data blocks (the
  "lower half") of all MassBank records for a compound,
  using the annotation data in the RMassBank options,
  spectrum info data from the \code{analyzedSpec}-type
  record and the peaks from the reanalyzed,
  multiplicity-filtered peak table. It calls
  \code{gatherSpectrum} for each child spectrum.
}
\details{
  The returned data blocks are in format \code{list(
  "AC\$MASS_SPECTROMETRY" = list('FRAGMENTATION_MODE' =
  'CID', ...), ...)} etc.
}
\note{
  Note that the global table \code{additionalPeaks} is also
  used as an additional source of peaks.
}
\examples{
\dontrun{
     myspectrum <- aggregatedRcSpecs$specComplete[[1]]
		massbankdata <- gatherCompound(myspectrum, refilteredRcSpecs)
		# Note: ac_lc and ac_ms are data blocks usually generated in gatherCompound and
		# passed on from there. The call below gives a relatively useless result :)
		ac_lc_dummy <- list()
		ac_ms_dummy <- list()
		justOneSpectrum <- gatherSpectrum(myspectrum, myspectrum$msmsdata[[2]],
			ac_ms_dummy, ac_lc_dummy, refilteredRcSpecs)
}
}
\author{
  Michael Stravs
}
\references{
  MassBank record format:
  \url{http://www.massbank.jp/manuals/MassBankRecord_en.pdf}
}
\seealso{
  \code{\link{mbWorkflow}}, \code{\link{compileRecord}}
}

