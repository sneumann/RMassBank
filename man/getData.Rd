% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpectrumMethods.R
\docType{methods}
\name{getData}
\alias{getData}
\alias{getData,RmbSpectrum2-method}
\title{Get data frame with all present peak data}
\usage{
\S4method{getData}{RmbSpectrum2}(s)
}
\arguments{
\item{s}{The \code{RmbSpectrum2} object to extract data from.}
}
\value{
A data frame with columns for every set slot.
}
\description{
Returns a data frame with columns for all non-empty slots in a \code{RmbSpectrum2} object. Note that \code{MSnbase::Spectrum} has
a method \code{as.data.frame}, however that one will return only mz, intensity. This function is kept separate to ensure downwards
compatibility since it returns more columns than MSnbase \code{as.data.frame}.
}
\author{
stravsmi
}
