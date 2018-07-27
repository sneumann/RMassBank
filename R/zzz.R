.onLoad <- function(libname, pkgname) {
  RMassBank.env <<- new.env()
  RMassBank.env$ReadAnnotation <- FALSE
  RMassBank.env$testnumber <- 1
  ## new variables
  RMassBank.env$verbose.output <- FALSE
  RMassBank.env$export.invalid <- FALSE
  RMassBank.env$export.molfiles <- TRUE
  
  mb <- list()
  attach(RMassBank.env)
}

utils::globalVariables(c("cpdID", "isotopes","mzCalc"))