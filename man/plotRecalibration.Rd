\name{plotRecalibration}
\alias{plotRecalibration}
\alias{plotRecalibration.direct}
\title{Plot the recalibration graph.}
\usage{
  plotRecalibration(w)

  plotRecalibration.direct(rcdata, rc, rc.ms1, title,
    mzrange)
}
\arguments{
  \item{w}{The workspace to plot the calibration graph
  from}

  \item{rcdata}{A data frame with columns \code{recalfield}
  and \code{mzFound}.}

  \item{rc}{Predictor for MS2 data}

  \item{rc.ms1}{Predictor for MS1 data}

  \item{title}{Prefix for the graph titles}

  \item{mzrange}{m/z value range for the graph}
}
\description{
  Plot the recalibration graph.
}
\author{
  Michele Stravs, Eawag <michael.stravs@eawag.ch>
}

