

# Accessor methods for found, complete, empty
# .checkSpectra <-
#' Check if a spectra set is found, complete, empty
#' 
#' Checks if a specific compound (\code{RmbSpectraSet}) was found with child spectra in the raw file (\code{found}), 
#' has a complete set of MS2 spectra with useful peaks (\code{complete}), or is empty (note: spectra are currently not ever
#' marked empty - empty should mean found, but no useful peaks at all. This is not yet currently tested.)
#' 
#' @param s The  (\code{RmbSpectraSet}) to check
#' @param property The property to check (\code{found}, \code{complete} or \code{empty})
#' @return \code{TRUE} or \code{FALSE}
#' 
#' @author stravsmi
#' @export
setGeneric("checkSpectra", function(s, property) standardGeneric("checkSpectra"))

# .checkSpectra.RmbSpectraSet <-
#' @describeIn checkSpectra 
#' @export
setMethod("checkSpectra", c("RmbSpectraSet", "character"), function(s, property)
		{
			#stopifnot(value=="logical", "For single spectraSet, only TRUE/FALSE output is supported.")
			fields <- c("found", "complete", "empty")
			if(!(property %in% fields)) stop("Only found, complete, empty properties are allowed")
			isTRUE(slot(s, property))
		})


# .selectSpectra <- 

#' Select a subset of spectra matching properties
#' 
#' From a list of \code{RmbSpectraSet}s, returns the spectra which match a criterion (found, complete, empty as in \code{\link{checkSpectra}}).
#' This can be returned either as a \code{TRUE/FALSE} vector, as a vector of indices for matching elements, as a vector of \code{RmbSpectraSet} objects
#' matching the conditions, or as a vector of \code{RmbSpectraSet} objects NOT matching the conditions (sic!).
#' 
#' @param s The \code{RmbSpectraSetList} or \code{msmsWorkspace} to select \code{RmbSpectraSet}s from. 
#' @param property The property to check (\code{found}, \code{complete} or \code{empty})
#' @param value \code{logical} if a \code{TRUE/FALSE} list should be returned; \code{index} if a vector of matching indices should be returned,
#' 			\code{object} if matching objects should be returned, \code{mismatch} if mismatching objects should be returned.
#' @return As described above.
#' 
#' @author stravsmi
#' @export
setGeneric("selectSpectra",  def = function(s, property, value="logical") standardGeneric("selectSpectra"),
	signature = c("s", "property"))

# .selectSpectra.RmbSpectraSetList <- 
#' @describeIn selectSpectra A method for selecting spectra from a spectra set list
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
#' @describeIn selectSpectra A method for selecting spectra from an msmsWorkspace
#' @export
setMethod("selectSpectra", c("msmsWorkspace", "character"), function(s, property, value="logical") 
			selectSpectra(s@spectra, property, value))

# .spectraCount <- 
#' Count MS2 spectra per compound
#' 
#' Counts the number of acquired spectra for a compound or multiple compounds
#' 
#' @param s  The object (\code{RmbSpectraSet}, \code{RmbSpectraSetList} or \code{msmsWorkspace}) to count the spectra in.
#' @return For \code{RmbSpectraSet} objects, a single number counting the spectra in that object. For \code{RmbSpectraSetList} or \code{msmsWorkspace}, a 
#' vector with spectra counts for all compounds (\code{RmbSpectraSet}s) in the object.
#' 
#' @author stravsmi
#' @export
setGeneric("spectraCount", function(s) standardGeneric("spectraCount"))

# .spectraCount.RmbSpectraSet <-
#' @describeIn spectraCount Counts the number of acquired spectra for an RmbSpectraSet
#' @export
setMethod("spectraCount", c("RmbSpectraSet"), function(s)
		{
			length(s@children)
		})

# .spectraCount.RmbSpectraSetList <- 
#' @describeIn spectraCount Counts the number of acquired spectra for an RmbSpectraSetList
#' @export
setMethod("spectraCount", c("RmbSpectraSetList"), function(s)
		{
			unlist(lapply(s, spectraCount))
		})


# .spectraCount.msmsWorkspace <-
#' @describeIn spectraCount Counts the number of acquired spectra for an msmsWorkSpace
#' @export
setMethod("spectraCount", c("msmsWorkspace"), function(s) 
			spectraCount(s@spectra))


#' @export
#' @describeIn selectPeaks A method to filter spectra to the specified peaks
setMethod("selectPeaks", c("RmbSpectraSetList"), function(o, ..., enclos = parent.frame(2))
		{
			for(n in seq_len(length(o)))
				o[[n]] <- selectPeaks(o[[n]], ..., enclos=enclos)
			return(o)
		})

#' @export 
setMethod("selectPeaks", c("RmbSpectraSet"), function(o, ..., enclos = parent.frame(2))
		{
			o@children <- selectPeaks(o@children, ..., enclos=enclos)
			return(o)
		})