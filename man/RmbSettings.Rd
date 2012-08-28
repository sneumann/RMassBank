\name{RmbSettings}
\alias{RmbSettings}
\title{RMassBank settings}
\description{
  Describes all settings for the RMassBank settings file.
}
\details{
  \itemize{ \item{\code{deprofile}}{Whether and how to
  deprofile input raw files. Leave the setting empty if
  your raw files are already in "centroid" mode. If your
  input files are in profile mode, you have the choice
  between algorithms \code{\link{deprofile}.spline,
  deprofile.fwhm, deprofile.localMax}; refer to the
  individual manpages for more information.}
  \item{\code{rtMargin, rtShift}}{The allowed retention
  time deviation relative to the values specified in your
  compound list (see \code{\link{loadList}}), and the
  systematic shift (due to the use of, e.g., pre-columns or
  other special equipment.} \item{\code{babeldir}}{
  Directory to OpenBabel. Required for creating molfiles
  for MassBank export.  If no OpenBabel directory is given,
  RMassBank will attempt to use the CACTUS webservice for
  SDF generation. It is strongly advised to install
  OpenBabel; the CACTUS structures have explicit hydrogen
  atoms.  The path should point to the directory where
  babel.exe (or the Linux "babel" equivalent) lies.  }
  \item{\code{use_version}}{Which MassBank record format to
  use; version 2 is strongly advised, version 1 is
  considered outdated and should be used only if for some
  reason you are running old servers and an upgrade is not
  feasible.} \item{\code{use_rean_peaks}}{Whether to
  include peaks from reanalysis (see
  \code{\link{reanalyzeFailpeaks}}) in the MassBank
  records. Boolean, TRUE or FALSE.  }
  \item{\code{annotations}}{ A list of constant annotations
  to use in the MassBank records. The entries
  \code{authors, copyright, license, instrument,
  instrument_type, compound_class} correspond to the
  MassBank entries \code{AUTHORS, COPYRIGHT, LICENSE,
  AC$INSTRUMENT, AC$INSTRUMENT_TYPE, CH$COMPOUND_CLASS}.
  The entry \code{confidence_comment} is added as
  \code{COMMENT: CONFIDENCE} entry.

  The entry \code{internal_id_fieldname} is used to name
  the MassBank entry which will keep a reference to the
  internal compound ID used in the workflow: for
  \code{internal_id_fieldname = MYID} and e.g. compound
  1234, an entry will be added to the MassBank record with
  \code{COMMENT: MYID 1234}. The internal fieldname should
  not be left empty!

  The entries \code{lc_gradient, lc_flow, lc_solvent_a,
  lc_solvent_b, lc_column} correspond to the MassBank
  entries \code{AC$CHROMATOGRAPHY: FLOW_GRADIENT,
  FLOW_RATE, SOLVENT A, SOLVENT B, COLUMN_NAME}.

  \code{ms_type, ionization} correspond to
  \code{AC$MASS_SPECTROMETRY: MS_TYPE, IONIZATION}.

  \code{entry_prefix} is the two-letter prefix used when
  building MassBank accession codes.

  Entries under \code{ms_dataprocessing} are added as
  \code{MS$DATA_PROCESSING:} entries, in addition to the
  default \code{WHOLE: RMassBank}.  }
  \item{\code{spectraList}}{This setting describes the
  experimental annotations for the single data-dependent
  scans. For every data-dependent scan event, a
  \code{spectraList} entry with \code{mode, ces, ce, res}
  denoting collision mode, collision energy in short and
  verbose notation, and FT resolution.}
  \item{\code{accessionNumberShifts}}{This denotes the
  starting points for accession numbers for different ion
  types. For example, \code{pH: 0, mH: 50} means that
  [M+H]+ spectra will start at \code{XX123401} (\code{XX}
  being the \code{entry_prefix} and \code{1234} the
  compound id) and [M-H]- will start at \code{XX123451}.}
  \item{\code{electronicNoise, electronicNoiseWidth}}{Known
  electronic noise peaks and the window to be used by
  \code{\link{cleanElnoise}}}
  \item{\code{recalibrateBy}}{\code{dppm} or \code{dmz} to
  recalibrate either by delta ppm or by delta mz.}
  \item{\code{recalibrateMS1}}{\code{common} or
  \code{separate} to recalibrate MS1 data points together
  or separately from MS2 data points.}
  \item{\code{recalibrator: MS1, MS2}}{The functions to use
  for recalibration of MS1 and MS2 data points.  Note that
  the \code{MS1} setting is only meaningful if
  \code{recalibrateMS1: separate}, otherwise the \code{MS2}
  setting is used for a common recalibration curve. See
  \code{\link{recalibrate.loess}} for details.  } }
}
\seealso{
  \code{\link{loadRmbSettings}}
}

