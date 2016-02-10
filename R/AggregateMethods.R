
#' Select matching/unmatching peaks from aggregate table
#' 
#' @param o Workspace or aggregate table from a workspace
#' @return Selects the peaks from the aggregate table which matched within filter criteria (\code{peaksMatched}) or didn't match
#' 		(\code{peaksUnmatched}).
#' 
#' @author stravsmi
#' @export
setGeneric("peaksMatched", function(o){ standardGeneric("peaksMatched")})

#' @export
#' @describeIn peaksMatched A method to retrieve the matched peaks from the "aggregated" slot (a data.frame object) in an msmsWorkSpace
setMethod("peaksMatched", c("data.frame"), function(o)
		{
			o[o$good,,drop=FALSE]
		})

#' @export
#' @describeIn peaksMatched A method to retrieve the matched peaks from an msmsWorkSpace
setMethod("peaksMatched", c("msmsWorkspace"), function(o) peaksMatched(o@aggregated))

#' Select matching/unmatching peaks from aggregate table
#' 
#' @param o Workspace or aggregate table from a workspace
#' @param cleaned Return only peaks which pass electronic noise filtering if \code{TRUE}.
#' @return Selects the peaks from the aggregate table which matched within filter criteria (\code{peaksMatched}) or didn't match
#' 		(\code{peaksUnmatched}).
#' 
#' @author stravsmi
#' @export
setGeneric("peaksUnmatched", function(o, cleaned=FALSE) standardGeneric("peaksUnmatched"))

#' @export
#' @describeIn peaksUnmatched A method to retrieve the unmatched peaks from the "aggregated" slot (a data.frame object) in an msmsWorkSpace
setMethod("peaksUnmatched", c("data.frame"), function(o, cleaned=FALSE)
		{
			if(cleaned)
			{
				if(!("noise" %in% colnames(o)))
					stop("Electronic noise cleaning hasn't been performed yet.")
				return(o[(!o$good) & (!o$noise),,drop=FALSE])
			}
			else
				o[!o$good,,drop=FALSE]
		})

#' @export
#' @describeIn peaksUnmatched A method to retrieve the unmatched peaks from an msmsWorkSpace
setMethod("peaksUnmatched", c("msmsWorkspace"), function(o, cleaned=FALSE) peaksUnmatched(o@aggregated, cleaned))



##' @export
##' @describeIn selectPeaks A method to retrieve the specified peaks from the "aggregated" slot (a data.frame object) in an msmsWorkSpace
#setMethod("selectPeaks", c("data.frame"), function(o, good=FALSE, bad=FALSE, cleaned=FALSE, best=FALSE)
#		{
#			if(!("noise" %in% colnames(o)))
#				o$noise <- FALSE
#			
#			o[
#					((good & o$good) |
#					 (bad & !o$good)) &
#				 (!cleaned | !o$noise) &
#				 (!best | (o$dppm == o$dppmBest))
#					,,drop=FALSE]			
#		})



#' @export
setMethod("selectPeaks", c("data.frame"), function(o, filter, ..., enclos=parent.frame(2))
    {
      if(missing(filter))
        return(o)
      
      f <- substitute(filter)
      o <- o[eval(f, o, enclos) & !is.na(eval(f,o,enclos)),,drop=FALSE]
    })

#' @export
#' @describeIn selectPeaks A method to retrieve the specified peaks from an msmsWorkSpace
setMethod("selectPeaks", c("msmsWorkspace"), 
    function(o, ..., enclos=parent.frame(2)) selectPeaks(o@aggregated, ..., enclos))


#' @export
#' @describeIn addProperty Add a new column to a data.frame
setMethod("addProperty", c("data.frame", "character", "character", "ANY"), function(o, name, type, value=NA)
		{
			o[,name] <- as(rep(value, nrow(o)), type)
			o
		})

