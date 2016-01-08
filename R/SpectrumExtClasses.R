# TODO: Add comment
# 
# Author: stravsmi
###############################################################################



#' @exportClass RmbSpectrum2Ext
.RmbSpectrum2Ext <- setClass("RmbSpectrum2Ext",
		representation = representation(
				properties="data.frame"
		),
		contains=c("RmbSpectrum2"),
		prototype = prototype(
				properties = data.frame(),
				new("Versioned", versions=c(classVersion("RmbSpectrum2"), RmbSpectrum2Ext = "0.1.0"))
		),
)


#' @describeIn addProperty Add a new column to the RmbSpectrum2 properties
#'
#' @export 
setMethod("addProperty", c("RmbSpectrum2", "character", "character", "ANY"), function(o, name, type, value=NA)
		{
			if(ncol(o@properties) == 0)
				o@properties <- data.frame(row.names = seq_len(o@peaksCount)) 
			o@properties[,name] <- as(rep(value, o@peaksCount), type)
			o
		})

#setGeneric("setData",	function(s, df, ...) standardGeneric("setData"))
setGeneric("property", function(o, property) standardGeneric("property"))
setGeneric("property<-", function(o, property, value, addNew, class) standardGeneric("property<-"))

#' @export
setMethod("property", c("RmbSpectrum2Ext", "character"), function(o, property)
		{
			if(property %in% colnames(o@properties))
				return(o@properties[,property])
			else
				# We can't use FALSE or NA, since it could be confused with a 1-length logical FALSE or 1-length ANY NA  
				return(NULL)
		})

#' @export
setMethod("property<-", c("RmbSpectrum2Ext", "character", "ANY", "logical", "character"), function(o, property, value, addNew = FALSE, class="")
		{
			if(class == "") class <- class(value)
			if(!(property %in% colnames(o@properties)) & !addNew)
			{
				warning("Trying to set inexistent property.")
				return(o)
			}
			else if(!(property %in% colnames(o@properties)))
				o <- addProperty(o, property, class)
			o@properties[,property] <- value
			return(o)
		})


#' @export 
setMethod("getData", c("RmbSpectrum2Ext"), function(s)
		{
			peaks <- s@peaksCount
			data <- callNextMethod()
			if(nrow(s@properties) == peaks)
				data <- cbind(data, s@properties)
			else if(nrow(s@properties) > 0)
				stop("Incorrect number of rows in properties frame.")
			return(data)
		})

#' @export
setMethod("setData", c("RmbSpectrum2Ext", "data.frame"), function(s, df, clean = TRUE)
		{
			# first set everything that RmbSpectrum2 setData does already
			s <- callNextMethod()
			# then find which columns can be set for properties and always clean (do no remove properties!)
			
			cols <- colnames(s@properties)
			cols.inDf <- which(cols %in% colnames(df))
			
			#newDf <- s@properties[rep(NA, s@peaksCount),,drop=FALSE]
			cols.df <- cols[cols.inDf]
			newDf <- df[,cols.df,drop=FALSE]
			
						
			cols.notinDf <- !(cols.inDf)
			cols.no <- cols[cols.notinDf]
			
			for(col in cols.no)
				newDf[,col] <- rep(as(NA, class(newDf[,col])), nrow(newDf))
			
			# reorder columns
			newDf <- newDf[,cols,drop=FALSE]
			
			s@properties <- newDf
			s
		})


#' @export
setMethod("initialize", "RmbSpectrum2Ext", function(.Object, ..., properties = .Object@properties) {
			## do work of initialization
			callNextMethod(.Object, ...,
					properties=properties)
		})



#' @export 
setMethod("normalize", c(object="RmbSpectrum2Ext"), function(object, ..., scale=999, precision=0, slot="intensity")
		{
			intensity <- .normalize.int(object, ...) * scale
			if(precision >= 0) intensity <- round(intensity, precision)
			if(slot != "intensity")
			{
				property(object, slot, TRUE, "numeric") <- intensity
			}
			else
				object@intensity <- intensity
			object
		})

## 
## #' @export
## RmbSpectrum2Ext <- function(s)
## {
##     sp <- new("RmbSpectrum2Ext")
## 
## }

setGeneric("spectrumUpgrade", function(o) standardGeneric("spectrumUpgrade"))

setMethod("spectrumUpgrade", "RmbSpectrum2List", function(o)
		{
			for(n in seq_len(length(o)))
				o[[n]] <- new("RmbSpectrum2Ext", o[[n]])
			return(o)
		})

setMethod("spectrumUpgrade", "RmbSpectraSet", function(o)
		{
			o@children <- spectrumUpgrade(o@children)
			return(o)
		})



setMethod("spectrumUpgrade", "msmsWorkspace", function(o)
		{
			o@spectra <- spectrumUpgrade(o@spectra)
			return(o)
		})
