# TODO: Add comment
# 
# Author: stravsmi
###############################################################################


setGeneric("mergePeaks",	function(peaks, ...) standardGeneric("mergePeaks"))
setGeneric("mergeSpectra", function(spectra, ...) standardGeneric("mergeSpectra"))

#' Merge peaks for spectra merging, FT shoulder elimination etc.
#' 
#' Note: ppm and abs are not cumulative!
#' @export
mergePeaks.df <- function(peaks, dppm, dabs, int)
{
	cutoff_int_limit <- int
	cutoff_mz_limit <- dabs
	cutoff_ppm_limit <- dppm
	# Order by intensity (descending)
	peaks_o <- peaks[order(peaks$intensity, decreasing=TRUE),,drop=FALSE]
	n <- 1
	# As long as there are peaks left AND the last peak is small enough (relative
	# to selected), move to the next peak
	while(n < nrow(peaks_o))
	{
		if(peaks_o[nrow(peaks_o),"intensity"] >= cutoff_int_limit *peaks_o[n,"intensity"])
			break
		# remove all peaks within cutoff_mz_limit (std. m/z = 0.5) which have intensity
		# of less than 5% relative to their "parent" peak
		#
		peaks_o <- peaks_o[ !(
							(
								((peaks_o$mz > peaks_o[n,"mz"] - cutoff_mz_limit)
									& (peaks_o$mz < peaks_o[n,"mz"] + cutoff_mz_limit))
								|  ((peaks_o$mz > peaks_o[n,"mz"] - ppm(peaks_o[n, "mz"], cutoff_ppm_limit, p=TRUE))
									& (peaks_o$mz < peaks_o[n,"mz"] + ppm(peaks_o[n, "mz"], cutoff_ppm_limit, p=TRUE)))
								)
							& (peaks_o$intensity < cutoff_int_limit * peaks_o[n,"intensity"])
							
							),,drop=FALSE]		 
		n <- n+1
	}
	return(peaks_o[order(peaks_o$mz),,drop=FALSE])
}

#' @export
setMethod("mergePeaks", "data.frame", function(peaks, ...)
		{
			mergePeaks.df(peaks, ...)
		})

#' @export
setMethod("mergePeaks", "matrix", function(peaks, ...)
		{
			mergePeaks.df(peaks, ...)
		})

#' @export
setMethod("mergePeaks", "RmbSpectrum2", function(peaks, ...)
		{
			df <- getData(peaks)
			df <- mergePeaks.df(df, ...)
			df <- setData(peaks, df)
			return(peaks)
		})

#' @export
setMethod("mergePeaks", "Spectrum", function(peaks, ...)
		{
			df <- as.data.frame(peaks)
			df <- mergePeaks.df(df, ...)
			peaks@mz <- df[,1]
			peaks@intensity <- df[,2]
			peaks@peaksCount <- nrow(df)
			return(peaks)
		})


#' @export
mergeSpectra.RmbSpectrum2List <- function(spectra)
{
	if(length(spectra) == 0)
		stop("Merging empty spectra lists doesn't work because of missing metadata!")
	df.all <- lapply(spectra, getData)
	df <- do.call(rbind, df.all)
	spectrum <- spectra[[1]]
	spectrum <- setData(spectrum, df)
	spectrum@collisionEnergy <- -1
	spectrum@merged <- length(spectra)
	return(spectrum)
}

#' @export
setMethod("mergeSpectra", "RmbSpectrum2List", function(spectra, ...)
			mergeSpectra.RmbSpectrum2List(spectra, ...))




