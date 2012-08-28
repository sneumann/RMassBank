\name{loadList}
\alias{loadList}
\alias{resetList}
\title{Load compound list for RMassBank}
\usage{
  loadList(path, listEnv=NULL)

  resetList()
}
\arguments{
  \item{path}{Path to the CSV list.}

  \item{listEnv}{The environment to load the list into. By
  default, the namelist is loaded into an environment
  internally in RMassBank.}
}
\value{
  No return value.
}
\description{
  Loads a CSV compound list with compound IDs
}
\details{
  The list is loaded into the variable
  \code{\var{compoundList}} in the environment
  \code{listEnv} (which defaults to the global environment)
  and used by the \code{findMz}, \code{findCAS}, ...
  functions.

  resetList() clears a currently loaded list.
}
\examples{
##
\dontrun{loadList("mylist.csv")}
}
\author{
  Michael Stravs
}
\seealso{
  \code{\link{findMz}}
}

