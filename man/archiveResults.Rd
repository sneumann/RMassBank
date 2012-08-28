\name{archiveResults}
\alias{archiveResults}
\title{Backup \code{msmsWorkflow} results}
\usage{
  archiveResults(w, fileName)
}
\arguments{
  \item{w}{The \code{msmsWorkspace} to be saved.}

  \item{fileName}{The filename to store the results under.}
}
\description{
  Writes the results from different \code{msmsWorkflow}
  steps to a file.
}
\examples{
# This doesn't really make a lot of sense,
		# it stores an empty workspace.
		w <- newMsmsWorkspace()
		archiveResults(w, "narcotics.RData")
}

