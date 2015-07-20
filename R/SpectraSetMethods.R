

# Accessor methods for found, complete, empty
# .checkSpectra <-
#' @export
setGeneric("checkSpectra", function(s, property) standardGeneric("checkSpectra"))

# .checkSpectra.RmbSpectraSet <-
#' @export
setMethod("checkSpectra", c("RmbSpectraSet", "character"), function(s, property)
		{
			#stopifnot(value=="logical", "For single spectraSet, only TRUE/FALSE output is supported.")
			fields <- c("found", "complete", "empty")
			if(!(property %in% fields)) stop("Only found, complete, empty properties are allowed")
			slot(s, property)
		})


# .selectSpectra <- 
#' @export
setGeneric("selectSpectra",  def = function(s, property, value="logical") standardGeneric("selectSpectra"),
	signature = c("s", "property"))

# .selectSpectra.RmbSpectraSetList <- 
#' @export
setMethod("selectSpectra", c("RmbSpectraSetList", "character"), function(s, property, value="logical")
		{
			matches <- unlist(lapply(s, function(s) checkSpectra(s, property)))
			if(value == "logical")
				return(matches)
			else if(value == "index")
				return(which(matches))	
			else if(value == "object")
				return(s[matches])
			else if(value == "mismatch")
				return(s[!matches])
		})

# .selectSpectra.msmsWorkspace <- 
#' @export
setMethod("selectSpectra", c("msmsWorkspace", "character"), function(s, property, value="logical") 
			selectSpectra(s@spectra, property, value))

# .spectraCount <- 
#' @export
setGeneric("spectraCount", function(s) standardGeneric("spectraCount"))

# .spectraCount.RmbSpectraSet <-
#' @export
setMethod("spectraCount", c("RmbSpectraSet"), function(s)
		{
			length(s@children)
		})

# .spectraCount.RmbSpectraSetList <- 
#' @export
setMethod("spectraCount", c("RmbSpectraSetList"), function(s)
		{
			unlist(lapply(s, spectraCount))
		})


# .spectraCount.msmsWorkspace <- 
#' @export
setMethod("spectraCount", c("msmsWorkspace"), function(s) 
			spectraCount(s@spectra))

