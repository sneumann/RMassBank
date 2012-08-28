\name{cleanElnoise}
\alias{cleanElnoise}
\title{Remove electronic noise}
\usage{
  cleanElnoise(peaks,
    noise=getOption("RMassBank")$electronicNoise, width =
    getOption("RMassBank")$electronicNoiseWidth)
}
\arguments{
  \item{peaks}{A data frame with peaks containing at least
  the columns \code{mzFound}, \code{dppm} and
  \code{dppmBest}.}

  \item{noise}{A numeric vector of known m/z of electronic
  noise peaks from the instrument Defaults to the entries
  in the RMassBank settings.}

  \item{width}{The window for the noise peak in m/z units.
  Defaults to the entries in the RMassBank settings.}
}
\value{
  Returns a dataframe where the rows matching electronic
  noise criteria are removed.  %% ...
}
\description{
  Removes known electronic noise peaks from a peak table
}
\examples{
# As used in the workflow:
\dontrun{
	    aggregatedRcSpecs$peaksUnmatchedC <-
		cleanElnoise(aggregatedRcSpecs$peaksUnmatched)
}
}
\author{
  Michael Stravs
}
\seealso{
  \code{\link{msmsWorkflow}}
}

