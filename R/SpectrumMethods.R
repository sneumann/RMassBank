# TODO: Add comment
# 
# Author: stravsmi
###############################################################################

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