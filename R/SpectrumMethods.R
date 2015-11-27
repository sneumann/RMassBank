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