# TODO: Add comment
# 
# Author: stravsmi
###############################################################################

#' Get data frame with all present peak data 
#' 
#' Returns a data frame with columns for all non-empty slots in a \code{RmbSpectrum2} object. Note that \code{MSnbase::Spectrum} has
#' a method \code{as.data.frame}, however that one will return only mz, intensity. This function is kept separate to ensure downwards
#' compatibility since it returns more columns than MSnbase \code{as.data.frame}.
#' @name getData
#' @aliases getData,RmbSpectrum2-method
#' 
#' @param s The \code{RmbSpectrum2} object to extract data from.
#' @return A data frame with columns for every set slot.
#' 
#' @author stravsmi
#' @docType methods
#' @export
setMethod("getData", c("RmbSpectrum2"), function(s)
		{
			peaks <- s@peaksCount
			cols <- c("mz", "intensity", "satellite", "low", "rawOK", "good", "mzCalc", "formula", "dbe", "formulaCount", "dppm", "dppmBest")
			cols.isFilled <- unlist(lapply(cols, function(col) length(slot(s, col)) == peaks))
			cols.filled <- cols[cols.isFilled]
			data <- lapply(cols.filled, function(col) slot(s, col))
			data$stringsAsFactors <- FALSE
			df <- do.call(data.frame,data)
			colnames(df) <- cols.filled
			df
		})




#' Set \code{RmbSpectrum2} data from data.frame
#' 
#' Sets all slots which are present as columns in the given dataframe. Optionally cleans the object, i.e. empties slots not defined in the data frame.
#' 
#' @name setData
#' @aliases setData,RmbSpectrum2,data.frame-method
#' 
#' @param s The \code{RmbSpectrum2} object to modify
#' @param df The data frame with new data
#' @param clean \code{TRUE} if slots which aren't present as columns in the data frame should be cleared.
#' @return The modified \code{RmbSpectrum2}.
#' 
#' @author stravsmi
#' @docType methods
#' @export
setMethod("setData", c("RmbSpectrum2", "data.frame"), function(s, df, clean = TRUE)
		{
			cols <- c("mz", "intensity", "satellite", "low", "rawOK", "good", "mzCalc", "formula", "dbe", "formulaCount", "dppm", "dppmBest")
			types <- c("mz" = "numeric", "intensity" = "numeric", "satellite" = "logical", "low" = "logical",
					"rawOK" = "logical", "good" = "logical", "mzCalc" = "numeric", "formula" = "character", 
					"dbe" = "numeric", "formulaCount" = "integer", "dppm" = "numeric", "dppmBest" = "numeric"
					)
			s@peaksCount <- as.integer(nrow(df))
			cols.inDf <- which(cols %in% colnames(df))
			cols.df <- cols[cols.inDf]
			for(col in cols.df)
			{
				slot(s, col) <- as(df[,col], types[[col]])
			}
			cols.notinDf <- !(cols.inDf)
			cols.no <- cols[cols.notinDf]
			if(clean)
			{
				for(col in cols.no)
				{
					slot(s, col) <- c()
				}
			}
			s
		})

setMethod("initialize", "RmbSpectrum2", function(.Object, ...,
				satellite = .Object@satellite,
				low = .Object@low,
				rawOK = .Object@rawOK,
				good = .Object@good,
				mzCalc = .Object@mzCalc,
				formula = .Object@formula,
				dbe = .Object@dbe,
				formulaCount = .Object@formulaCount,
				dppm = .Object@dppm,
				dppmBest = .Object@dppmBest,
				ok = .Object@ok,
				info = .Object@info
				) {
			## do work of initialization
			callNextMethod(.Object, ..., 
					satellite = satellite,
					low = low,
					rawOK = rawOK,
					good = good,
					mzCalc = mzCalc,
					formula = formula,
					dbe = dbe,
					formulaCount = formulaCount,
					dppm = dppm,
					dppmBest = dppmBest,
					ok = ok,
					info = info)
		})



#' @export
#' @describeIn selectPeaks A method to filter spectra to the specified peaks
setMethod("selectPeaks", c("RmbSpectrum2"), function(o, filter, ...)
		{
			if(missing(filter))
				return(o)
			df <- getData(o)
			f <- substitute(filter)
			df <- df[eval(f, df),,drop=FALSE]
			o <- setData(o, df)
			o
		})


setMethod("selectPeaks", c("Spectrum"), function(o, filter, ...)
		{
			if(missing(filter))
				return(o)
			df <- as.data.frame(o)
			f <- substitute(filter)
			df <- df[eval(f, df),,drop=FALSE]
			o@mz <- df[,1]
			o@intensity <- df[,2]
			o
		})


#' @export
#' @describeIn selectPeaks A method to filter spectra to the specified peaks
setMethod("selectPeaks", c("RmbSpectrum2List"), function(o, ...)
		{
			s <- lapply(o, function(s) selectPeaks(s, ...))
			for(n in seq_len(length(o)))
				o[[n]] <- s[[n]]
			return(o)
		})

#' @export 
setMethod("normalize", c(object="Spectrum"), function(object, ..., scale=999, precision=0)
		{
			intensity <- .normalize.int(object, ...) * scale
			if(precision >= 0) intensity <- round(intensity, precision)
			object@intensity <- intensity
			object
		})

.normalize.int <- function(object, ...)
{
	spec <- selectPeaks(object, ...)
	maxint <- max(spec@intensity)
	intensity <- object@intensity / maxint
	return(intensity)
}

#' @export
setMethod("normalize", c("RmbSpectrum2List"), function(object, ...)
		{
			s <- lapply(object, function(s) normalize(s, ...))
			for(n in seq_len(length(object)))
				object[[n]] <- s[[n]]
			return(object)
		})



setMethod("+", c("Spectrum", "numeric"), function(e1, e2) 
		{
			e1@mz <- e1@mz + e2
			return(e1)
		}) 



setMethod("-", c("Spectrum", "numeric"), function(e1, e2) 
		{
			e1@mz <- e1@mz - e2
			return(e1)
		}) 

setMethod("+", c("RmbSpectraSet", "ANY"), function(e1, e2)
		{
			e1@parent <- e1@parent + e2
			for(n in seq_len(length(e1@children)))
				e1@children[[n]] <- e1@children[[n]] + e2
			e1
		})

setMethod("-", c("RmbSpectraSet", "ANY"), function(e1, e2)
		{
			e1@parent <- e1@parent - e2
			for(n in seq_len(length(e1@children)))
				e1@children[[n]] <- e1@children[[n]] - e2
			e1
		})


setMethod("+", c("RmbSpectrum2List", "ANY"), function(e1, e2)
		{
			for(n in seq_len(length(e1)))
				e1[[n]] <- e1[[n]] + e2
			e1
		})


setMethod("-", c("RmbSpectrum2List", "ANY"), function(e1, e2)
		{
			for(n in seq_len(length(e1)))
				e1[[n]] <- e1[[n]] - e2
			e1
		})
