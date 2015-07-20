
#' @export
setGeneric("peaksMatched", function(o) standardGeneric("peaksMatched"))

#' @export
setMethod("peaksMatched", c("data.frame"), function(o)
		{
			o[o$good,,drop=FALSE]
		})

#' @export
setMethod("peaksMatched", c("msmsWorkspace"), function(o) peaksMatched(o@aggregated))

#' @export
setGeneric("peaksUnmatched", function(o) standardGeneric("peaksUnmatched"))

#' @export
setMethod("peaksUnmatched", c("data.frame"), function(o)
		{
			o[!o$good,,drop=FALSE]
		})

#' @export
setMethod("peaksUnmatched", c("msmsWorkspace"), function(o) peaksUnmatched(o@aggregated))