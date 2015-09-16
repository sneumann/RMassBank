.onLoad <- function(libname, pkgname) {
  RMassBank.env <<- new.env()
  RMassBank.env$ReadAnnotation <- FALSE
  RMassBank.env$testnumber <- 1
  mb <- list()
  attach(RMassBank.env)
}

utils::globalVariables(c("cpdID", "isotopes","mzCalc"))