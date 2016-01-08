.onLoad <- function(libname, pkgname) {
  RMassBank.env <<- new.env()
  RMassBank.env$ReadAnnotation <- FALSE
  RMassBank.env$testnumber <- 1
  mb <- list()
  attach(RMassBank.env)
}

utils::globalVariables(c("cpdID", "isotopes","mzCalc"))


# Overwrite the Versioned initialize function!

setMethod("initialize", "Versioned", 
		
		function (.Object, ...) 
		{
			.local <- function (.Object, ..., versions = list()) 
			{
				.Object <- callNextMethod(.Object, ...)
				classVersion(.Object)[names(versions)] <- versions
				.Object
			}
			.local(.Object, ...)
		})
