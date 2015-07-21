
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
setGeneric("peaksUnmatched", function(o, cleaned=FALSE) standardGeneric("peaksUnmatched"))

#' @export
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
setMethod("peaksUnmatched", c("msmsWorkspace"), function(o, cleaned=FALSE) peaksUnmatched(o@aggregated, cleaned))

#' @export 
setGeneric("selectPeaks", function(o, ...) standardGeneric("selectPeaks"))

#' @export
setMethod("selectPeaks", c("data.frame"), function(o, good=FALSE, bad=FALSE, cleaned=FALSE, best=FALSE)
		{
			if(!("noise" %in% colnames(o)))
				o$noise <- FALSE
			
			o[
					((good & o$good) |
					 (bad & !o$good)) &
				 (!cleaned | !o$noise) &
				 (!best | (o$dppm == o$dppmBest))
					,,drop=FALSE]			
		})

#' @export
setMethod("selectPeaks", c("msmsWorkspace"), function(o, ...) selectPeaks(o@aggregated, ...))

#' @export
setGeneric("addProperty", function(o, name, type, value=NA) standardGeneric("addProperty"))

setMethod("addProperty", c("data.frame", "character", "character", "ANY"), function(o, name, type, value=NA)
		{
			o[,name] <- as(rep(value, nrow(o)), type)
			o
		})



## #' @export 
## setGeneric("selectPeaks<-", function(o, ..., createColumns = TRUE, value) standardGeneric("selectPeaks<-"))
## 
## setMethod("selectPeaks<-", c("data.frame"), function(o, good=FALSE, bad=FALSE, cleaned=FALSE, best=FALSE, createColumns=TRUE, value)
##         {
##             if(!("noise" %in% colnames(o)))
##                 o$noise <- FALSE
## 
##             o[
##                     ((good & o$good) |
##                                 (bad & !o$good)) &
##                             (!cleaned | !o$noise) &
##                             (!best | (o$dppm == o$dppmBest))
##                     ,,drop=FALSE]			
##         })
## 
## setMethod("selectPeaks<-", c("msmsWorkspace"), function(o, ...,createColumns = TRUE, value) `selectPeaks<-`(o@aggregated, ..., createColumns, value))
## 


