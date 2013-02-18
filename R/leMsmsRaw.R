#' @import mzR
#' @importClassesFrom mzR
#' @importMethodsFrom mzR
# library(mzR)
NULL

#' Extract MS/MS spectra for specified precursor
#' 
#' Extracts MS/MS spectra from LC-MS raw data for a specified precursor, specified
#' either via the RMassBank compound list (see \code{\link{loadList}}) or via a mass.
#' 
#' Different versions of the function get the data from different sources.
#' 
#' @usage findMsMsHR(fileName, cpdID, mode="pH",confirmMode =0, useRtLimit = TRUE, dppm=10)
#' 
#' 		findMsMsHR.mass(msRaw, mz, limit.coarse, limit.fine, rtLimits = NA, maxCount = NA,
#' 		headerCache = NA)
#' 
#' 		findMsMsHR.direct(msRaw, cpdID, mode = "pH", confirmMode = 0,
#'  	useRtLimit = TRUE, dppm=10, limit.coarse=0.5)
#' 
#' @aliases findMsMsHR.mass findMsMsHR.direct findMsMsHR
#' @param fileName The file to open and search the MS2 spectrum in.
#' @param msRaw The opened raw file (mzR file handle) to search the MS2 spectrum in.
#' @param cpdID The compound ID in the compound list (see \code{\link{loadList}})
#' 			to use for formula lookup.
#' @param mz The mass to use for spectrum search.
#' @param dppm The limit in ppm to use for fine limit (see below) calculation.
#' @param limit.coarse The coarse limit to use for locating potential MS2 scans:
#'			this tolerance is used when finding scans with a suitable precursor
#' 			ion value.  
#' @param limit.fine The fine limit to use for locating MS2 scans: this tolerance
#' 			is used when locating an appropriate analyte peak in the MS1 precursor
#' 			spectrum.
#' @param mode The processing mode (determines which ion/adduct is searched):
#' 			\code{"pH", "pNa", "pM", "mH", "mM", "mFA"} for different ions 
#' 			([M+H]+, [M+Na]+, [M]+, [M-H]-, [M]-, [M+FA]-). 
#' @param confirmMode Whether to use the highest-intensity precursor (=0), second-
#' 			highest (=1), third-highest (=2)...
#' @param useRtLimit Whether to respect retention time limits from the compound list.
#' @param rtLimits \code{c(min, max)}: Minimum and maximum retention time to use
#' 			when locating the MS2 scans. 
#' @param headerCache If present, the complete \code{mzR::header(msRaw)}. Passing
#' 			this value is useful if spectra for multiple compounds should be 
#' 			extracted from the same mzML file, since it avoids getting the data
#' 			freshly from \code{msRaw} for every compound.
#' @param maxCount The maximal number of spectra groups to return. One spectra group
#' 			consists of all data-dependent scans from the same precursor whose precursor
#' 			mass matches the specified search mass.
#' @return	For \code{findMsMsHR} and \code{findMsMsHR.direct}: A "spectrum set", a list with items:
#' 			\item{foundOK}{\code{TRUE} if a spectrum was found, \code{FALSE} otherwise.
#' 				Note: if \code{FALSE}, all other values can be missing!}
#' 			\item{parentScan}{The scan number of the precursor scan.}
#' 			\item{parentHeader}{The header row of the parent scan, as returned by 
#' 				\code{mzR::header}.}
#' 			\item{childScans}{The scan numbers of the data-dependent MS2 scans.}
#' 			\item{childHeaders}{The header rows of the MS2 scan, as returned by
#' 				\code{mzR::header}.}
#' 			\item{parentPeak}{The MS1 precursor spectrum as a 2-column matrix}
#' 			\item{peaks}{A list of  2-column \code{mz, int} matrices of the MS2 scans.}
#' 			For \code{findMsMsHR.mass}: a list of "spectrum sets" as defined above, sorted
#' 			by decreasing precursor intensity.
#' 
#' @examples \dontrun{
#' 			loadList("mycompoundlist.csv")
#' 			# if Atrazine has compound ID 1:
#' 			msms_atrazine <- findMsMsHR("Atrazine_0001_pos.mzML", 1, "pH")
#' 			# Or alternatively:
#' 			msRaw <- openMSfile("Atrazine_0001_pos.mzML")
#' 			msms_atrazine <- findMsMsHR.direct(msRaw, 1, "pH")
#' 			# Or directly by mass (this will return a list of spectra sets):
#' 			mz <- findMz(1)$mzCenter
#' 			msms_atrazine_all <- findMsMsHR.mass(msRaw, mz, 1, ppm(msRaw, 10, p=TRUE))
#' 			msms_atrazine <- msms_atrazine_all[[1]]
#' }
#' @author Michael A. Stravs, Eawag <michael.stravs@@eawag.ch>
#' @seealso findEIC
#' @export
findMsMsHR <- function(fileName, cpdID, mode="pH",confirmMode =0, useRtLimit = TRUE, dppm=10)
{
	
	# access data directly for finding the MS/MS data. This is done using
	# mzR.
	msRaw <- openMSfile(fileName)
	ret <- findMsMsHR.direct(msRaw, cpdID, mode, confirmMode, useRtLimit, dppm)
	mzR::close(msRaw)
	return(ret)
}

#' @export
findMsMsHRperxcms <- function(fileName) {
	
	splitfn <- strsplit(fileName,'_')
    splitsfn <- splitfn[[1]]
    cpdID <- as.numeric(splitsfn[[length(splitsfn)-1]])
	
	parentMass <- findMass(cpdID) + 1
	RT <- findRt(cpdID)$RT * 60
	mzabs <- 0.1 ## TODO: als Parameter
	
	getRT <- function(xa) {
		rt <- sapply(xa@pspectra, function(x) {median(peaks(xa@xcmsSet)[x, "rt"])})
	}
	##
	## MSMS
	##
	xrmsms <- xcmsRaw(fileName, includeMSn=TRUE)

	## Where is the wanted isolation ?
	precursorrange <- range(which(xrmsms@msnPrecursorMz == parentMass)) ## TODO: add ppm one day

	## Fake MS1 from MSn scans
	xrmsmsAsMs <- msn2xcms(xrmsms)

	## Fake s simplistic xcmsSet
	xsmsms <-  xcmsSet (files=fileName,
					method="MS1")

	peaks(xsmsms) <- findPeaks(xrmsmsAsMs, method="centWave", peakwidth=c(5,12),
							prefilter=c(0,0), ppm=25, snthr=2,
							scanrange=precursorrange) ## TODO: parameters from findMsMsHRperxcms

	## Get pspec 
	pl <- peaks(xsmsms)[,c("mz", "rt")]
	candidates <- which( pl[,"mz"] < parentMass+mzabs & pl[,"mz"] > parentMass-mzabs
						& pl[,"rt"] < RT * 1.1 & pl[,"rt"] > RT * 0.9 )

	anmsms <- xsAnnotate(xsmsms)
	anmsms <- groupFWHM(anmsms)

	## Now find the pspec for Chelidonine
	psp <- which(sapply(anmsms@pspectra, function(x) {candidates %in% x}))
	
	## Alternative: Spectrum closest to MS1
	##psp <- which.min(getRT(anmsms) - actualRT)
	
	return(getpspectra(anmsms, psp))
}

#' @export
findMsMsHRperhand <- function(fileName) {
	
}

#' @export
findMsMsHR.mass <- function(msRaw, mz, limit.coarse, limit.fine, rtLimits = NA, maxCount = NA,
		headerCache = NA)
{
	eic <- findEIC(msRaw, mz, limit.fine, rtLimits)
	#	if(!is.na(rtLimits))
	#	{  
	#		eic <- subset(eic, rt >= rtLimits[[1]] & rt <= rtLimits[[2]])
	#	}
	if(!is.na(headerCache))
		headerData <- headerCache
	else
		headerData <- as.data.frame(header(msRaw))
	
	# Find MS2 spectra with precursors which are in the allowed 
	# scan filter (coarse limit) range
	findValidPrecursors <- headerData[
			(headerData$precursorMZ > mz - limit.coarse) &
			(headerData$precursorMZ < mz + limit.coarse),]
	# Find the precursors for the found spectra
	validPrecursors <- unique(findValidPrecursors$precursorScanNum)
	# check whether the precursors are real: must be within fine limits!
	# previously even "bad" precursors were taken. e.g. 1-benzylpiperazine
	which_OK <- lapply(validPrecursors, function(pscan)
			{
				pplist <- as.data.frame(
						mzR::peaks(msRaw, which(headerData$acquisitionNum == pscan)))
				colnames(pplist) <- c("mz","int")
				pplist <- pplist[(pplist$mz >= mz -limit.fine)
								& (pplist$mz <= mz + limit.fine),]
				if(nrow(pplist) > 0)
					return(TRUE)
				return(FALSE)
			})
	validPrecursors <- validPrecursors[which(which_OK==TRUE)]
	# Crop the "EIC" to the valid precursor scans
	eic <- eic[eic$scan %in% validPrecursors,]
	# Order by intensity, descending
	eic <- eic[order(eic$intensity, decreasing=TRUE),]
	if(nrow(eic) == 0)
		return(list(list(foundOK = FALSE)))
	if(!is.na(maxCount))
	{
		spectraCount <- min(maxCount, nrow(eic))
		eic <- eic[1:spectraCount,]
	}
	# Construct all spectra groups in decreasing intensity order
	spectra <- lapply(eic$scan, function(masterScan)
			{
				masterHeader <- headerData[headerData$acquisitionNum == masterScan,]
				childHeaders <- headerData[(headerData$precursorScanNum == masterScan) 
					& (headerData$precursorMZ > mz - limit.coarse) 
					& (headerData$precursorMZ < mz + limit.coarse) ,]
				childScans <- childHeaders$acquisitionNum
				
				msPeaks <- mzR::peaks(msRaw, masterHeader$seqNum)
				# if deprofile option is set: run deprofiling
				deprofile.setting <- getOption("RMassBank")$deprofile
				if(!is.na(deprofile.setting))
					msPeaks <- deprofile.scan(
							msPeaks, method = deprofile.setting, noise = NA, colnames = FALSE
							)
				colnames(msPeaks) <- c("mz","int")
				msmsPeaks <- lapply(childHeaders$seqNum, function(scan)
						{
							pks <- mzR::peaks(msRaw, scan)
							if(!is.na(deprofile.setting))
							{								
								pks <- deprofile.scan(
										pks, method = deprofile.setting, noise = NA, colnames = FALSE
								)
							}
							colnames(pks) <- c("mz","int")
							return(pks)
						}
				)
				return(list(
								foundOK = TRUE,
								parentScan = masterScan,
								parentHeader = masterHeader,
								childScans = childScans,
								childHeaders= childHeaders,
								parentPeak=msPeaks,
								peaks=msmsPeaks
						#xset=xset#,
						#msRaw=msRaw
						))
			})
	names(spectra) <- eic$acquisitionNum
	return(spectra)
}

#' @export
findMsMsHR.direct <- function(msRaw, cpdID, mode = "pH", confirmMode = 0, useRtLimit = TRUE, dppm=10, limit.coarse=0.5)
{
  # for finding the peak RT: use the gauss-fitted centwave peak
  # (centroid data converted with TOPP is necessary. save as
  # mzData, since this is correctly read :P)
  #xset <- xcmsSet(fileName, method="centWave",ppm=5, fitgauss=TRUE)

  # find cpd m/z
  mzLimits <- findMz(cpdID, mode)
  mz <- mzLimits$mzCenter
  limit.fine <- ppm(mz, dppm, p=TRUE)
  if(!useRtLimit)
	  rtLimits <- NA
  else
  {
	  rtMargin <- getOption("RMassBank")$rtMargin
	  dbRt <- findRt(cpdID)
	  rtLimits <- c(dbRt$RT - rtMargin, dbRt$RT + rtMargin) * 60
  }
  spectra <- findMsMsHR.mass(msRaw, mz, limit.coarse, limit.fine, rtLimits, confirmMode + 1)
  spectra[[confirmMode + 1]]$mz <- mzLimits
  return(spectra[[confirmMode + 1]])
}


# Finds the EIC for a mass trace with a window of x ppm.
# (For ppm = 10, this is +5 / -5 ppm from the non-recalibrated mz.)
#' Extract EICs 
#' 
#' Extract EICs from raw data for a determined mass window.
#' 
#' @param msRaw The mzR file handle 
#' @param mz The mass or mass range to extract the EIC for: either a single mass
#' 			(with the range specified by \code{limit} below) or a mass range
#' 			in the form of \code{c(min, max)}. 
#' @param limit If a single mass was given for \code{mz}: the mass window to extract.
#' 			A limit of 0.001 means that the EIC will be returned for \code{[mz - 0.001, mz + 0.001]}.
#' @param rtLimit If given, the retention time limits in form \code{c(rtmin, rtmax)} in seconds.
#' @return A \code{[rt, intensity, scan]} matrix (\code{scan} being the scan number.) 
#' @author Michael A. Stravs, Eawag <michael.stravs@@eawag.ch>
#' @seealso findMsMsHR
#' @export
findEIC <- function(msRaw, mz, limit = NULL, rtLimit = NA)
{
	# calculate mz upper and lower limits for "integration"
	if(all(c("mzMin", "mzMax") %in% names(mz)))
		mzlimits <- c(mz$mzMin, mz$mzMax)
	else
		mzlimits <- c(mz - limit, mz + limit)
	# Find peaklists for all MS1 scans
	headerData <- as.data.frame(header(msRaw))
	# If RT limit is already given, retrieve only candidates in the first place,
	# since this makes everything much faster.
	if(all(!is.na(rtLimit)))
		headerMS1 <- headerData[
				(headerData$msLevel == 1) & (headerData$retentionTime >= rtLimit[[1]])
						& (headerData$retentionTime <= rtLimit[[2]])
				,]
	else
		headerMS1 <- headerData[headerData$msLevel == 1,]
	pks <- mzR::peaks(msRaw, headerMS1$seqNum)
	# Sum intensities in the given mass window for each scan
	pks_t <- unlist(lapply(pks, function(peaktable)
						sum(peaktable[which((peaktable[,1] >= mzlimits[[1]]) & (peaktable[,1] <= mzlimits[[2]])) ,2])))
	rt <- headerMS1$retentionTime
	scan <- headerMS1$acquisitionNum
	return(data.frame(rt = rt, intensity=pks_t, scan=scan))
}


if (FALSE) {
library(RMassBank)
loadList("/home/sneumann/CASMI/example/Chelidonine.csv")


}
