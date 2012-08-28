\name{findMz.formula}
\alias{findMz.formula}
\title{Find the exact mass +/- a given margin for a given formula or its ions and adducts.}
\usage{
  findMz.formula(formula, mode = "pH", ppm = 10,
    deltaMz = 0)
}
\arguments{
  \item{formula}{The molecular formula in text or list
  format (see \code{\link{formulastring.to.list}}}

  \item{mode}{\code{"pH", "pNa", "pM", "mH", "mM", "mFA"}
  for different ions ([M+H]+, [M+Na]+, [M]+, [M-H]-, [M]-,
  [M+FA]-). "" for the uncharged molecule.}

  \item{ppm}{The ppm margin to add/subtract}

  \item{deltaMz}{The absolute mass to add/subtract.
  Cumulative with \code{ppm}}
}
\value{
  A \code{list(mzMin=, mzCenter=, mzMax=)} with the masses.
}
\description{
  Find the exact mass +/- a given margin for a given
  formula or its ions and adducts.
}
\examples{
findMz.formula("C6H6")
}
\author{
  Michael A. Stravs, Eawag <michael.stravs@eawag.ch>
}
\seealso{
  \code{\link{findMz}}
}

