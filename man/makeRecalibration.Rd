\name{makeRecalibration}
\alias{makeRecalibration}
\alias{recalibrateSingleSpec}
\alias{recalibrateSpectra}
\title{Recalibrate MS/MS spectra}
\usage{
  makeRecalibration(spec, mode)

  recalibrateSpectra(mode, rawspec = NULL, rc = NULL,
    rc.ms1 = NULL, w = NULL)

  recalibrateSingleSpec(spectrum, rc)
}
\arguments{
  \item{spec}{For \code{recalibrateSpectra}: a list of
  \code{aggregatedSpecs} type (i.e. as returned by
  \code{aggregateSpectra}).}

  \item{spectrum}{For \code{recalibrateSingleSpec}: a
  matrix with columns \code{mz, int} to be recalibrated.}

  \item{mode}{\code{"pH", "pNa", "pM", "mH", "mM", "mFA"}
  for different ions ([M+H]+, [M+Na]+, [M]+, [M-H]-, [M]-,
  [M+FA]-).}

  \item{rawspec}{For \code{recalibrateSpectra}:a
  \code{list} of \code{specs}-type object, i.e. as returned
  by the \code{\link{findMsMsHR}} function family.  If
  empty, no spectra are recalibrated, but the recalibration
  curve is returned.}

  \item{rc,rc.ms1}{The recalibration curves to be used in
  the recalibration.}

  \item{w}{The \code{msmsWorkspace} to write the
  calibration to or to get the calibration from.}
}
\value{
  \code{makeRecalibration}: a \code{list(rc, rc.ms1)} with
  recalibration curves for the MS2 and MS1 spectra.

  \code{recalibrateSpectra}: if \code{rawspec} is not
  \code{NULL}, returns the recalibrated spectra in the same
  structure as the input spectra. Each spectrum matrix has
  an additional column \code{mzRecal} with the recalibrated
  mass.

  \code{recalibrateSingleSpec}: a matrix with the single
  recalibrated spectrum.  Column \code{mzRecal} contains
  the recalibrated value.
}
\description{
  Recalibrates MS/MS spectra by building a recalibration
  curve of the assigned putative fragments of all spectra
  in \code{aggregatedSpecs} (measured mass vs. mass of
  putative associated fragment) and additionally the parent
  ion peaks.
}
\details{
  Note that the actually used recalibration functions are
  governed by the general MassBank settings (see
  \code{\link{recalibrate}}).

  If a set of acquired LC-MS runs contains spectra for two
  different ion types (e.g. [M+H]+ and [M+Na]+) which
  should both be processed by RMassBank, it is necessary to
  do this in two separate runs. Since it is likely that one
  ion type will be the vast majority of spectra (e.g. most
  in [M+H]+ mode), and only few spectra will be present for
  other specific adducts (e.g. only few [M+Na]+ spectra),
  it is possible that too few spectra are present to build
  a good recalibration curve using only e.g. the [M+Na]+
  ions. Therefore we recommend, for one set of LC/MS runs,
  to build the recalibration curve for one ion type
  (\code{msmsWorkflow(mode="pH", steps=c(1:8),
  newRecalibration=TRUE)}) and reuse the same curve for
  processing different ion types
  (\code{msmsWorkflow(mode="pNa", steps=c(1:8),
  newRecalibration=FALSE)}). This also ensures a consistent
  recalibration across all spectra of the same batch.
}
\examples{
\dontrun{
			rcCurve <- recalibrateSpectra(aggregatedSpecs, "pH")
			recalibratedSpecs <- recalibrateSpectra(aggregatedSpecs, "pH", specs, w=myWorkspace)
			recalibratedSpecs <- recalibrateSpectra(aggregatedSpecs, "pH", specs,
				rcCurve$rc, rcCurve$rc.ms1)
			s <- matrix(c(100,150,200,88.8887,95.0005,222.2223), ncol=2)
			colnames(s) <- c("mz", "int")
			recalS <- recalibrateSingleSpec(s, rcCurve$rc)
			}
}
\author{
  Michael Stravs, Eawag <michael.stravs@eawag.ch>
}

