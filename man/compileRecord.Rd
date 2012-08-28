\name{compileRecord}
\alias{compileRecord}
\title{Compile MassBank records}
\usage{
  compileRecord(spec, mbdata, refiltered, additionalPeaks =
    NULL)
}
\arguments{
  \item{spec}{A spectra block for a compound, as returned
  from \code{\link{analyzeMsMs}}. Note that \bold{peaks are
  not read from this object anymore}: Peaks come from the
  \code{refiltered} dataframe (and from the global
  \code{additionalPeaks} dataframe; cf.
  \code{\link{addPeaks}} for usage information.)}

  \item{mbdata}{The information data block for the record
  header, as stored in \code{mbdata_relisted} after loading
  an infolist.}

  \item{refiltered}{A list with at least the member
  \code{peaksOK}, and if peaks from reanalysis should be
  used, also \code{peaksReanOK}. \code{peaksOK} must be a
  dataframe with at least the, containing at least the
  columns \code{cpdID, scan, mzFound, formula, int, dppm}.
  If reanalyzed peaks are used, the column setup of
  \code{peaksReanOK} must be such as returned from
  \code{\link{filterMultiplicity}}.}

  \item{additionalPeaks}{If present, a table with
  additional peaks to add into the spectra.  As loaded with
  \code{\link{addPeaks}}.}
}
\value{
  Returns a MassBank record in list format: e.g.
  \code{list("ACCESSION" = "XX123456", "RECORD_TITLE" =
  "Cubane", ..., "CH\$LINK" = list( "CAS" = "12-345-6",
  "CHEMSPIDER" = 1111, ...))}
}
\description{
  Takes a spectra block for a compound, as returned from
  \code{\link{analyzeMsMs}}, and an aggregated cleaned peak
  table, together with a MassBank information block, as
  stored in the infolists and loaded via
  \code{\link{loadInfolist}}/\code{\link{readMbdata}} and
  processes them to a MassBank record
}
\details{
  \code{compileRecord} calls \code{\link{gatherCompound}}
  to create blocks of spectrum data, and finally fills in
  the record title and accession number, renames the
  "internal ID" comment field and removes dummy fields.
}
\examples{
#
\dontrun{myspec <- aggregatedRcSpecs\\$specFound[[1]]}
# after having loaded an infolist:
\dontrun{mbdata <- mbdata_relisted[[which(mbdata_archive\\$id == as.numeric(myspec\\$id))]]}
\dontrun{compiled <- compileRecord(myspec, mbdata, reanalyzedRcSpecs)}
}
\author{
  Michael Stravs
}
\references{
  MassBank record format:
  \url{http://www.massbank.jp/manuals/MassBankRecord_en.pdf}
}
\seealso{
  \code{\link{mbWorkflow}}, \code{\link{addPeaks}},
  \code{\link{gatherCompound}}, \code{\link{toMassbank}}
}

