#' 
#' Extracts and processes spectra from a specified file list, according to 
#' loaded options and given parameters.
#' 
#' The filenames of the raw LC-MS runs are read from the array \code{files} 
#' in the global enviroment.
#' See the vignette \code{vignette("RMassBank")} for further details about the
#' workflow.
#' 
#' @usage msmsRead(w, filetable = NULL, files = NULL, cpdids = NULL, 
#'					readMethod, mode, confirmMode = FALSE, useRtLimit = TRUE, 
#'					Args, settings = getOption("RMassBank"), progressbar = "progressBarHook")
#' @param w A \code{msmsWorkspace} to work with.
#' @param filetable The path to a .csv-file that contains the columns "files" and "cpdid" supplying
#' 			the relationships between files and compound IDs. Either this or "files" need
#'			to be specified. 
#' @param files A vector or list containing the filenames of the files that are to be read as spectra. 
#'				For the IDs to be inferred from the filenames alone, there need to be exactly 2 underscores.
#' @param cpdids A vector or list containing the compound IDs of the files that are to be read as spectra.
#'				The ordering of this and \code{files} implicitly assigns each ID to the corresponding file.
#'				If this is supplied, then the IDs implicitly named in the filenames are ignored.
#' @param readMethod Several methods are available to get peak lists from the files.
#'        Currently supported are "mzR", "xcms", "MassBank" and "peaklist".
#'        The first two read MS/MS raw data, and differ in the strategy 
#'        used to extract peaks. MassBank will read existing records, 
#'        so that e.g. a recalibration can be performed, and "peaklist" 
#'        just requires a CSV with two columns and the column header "mz", "int".
#' @param mode \code{"pH", "pNa", "pM", "mH", "mM", "mFA"} for different ions 
#' 			([M+H]+, [M+Na]+, [M]+, [M-H]-, [M]-, [M+FA]-).
#' @param confirmMode Defaults to false (use most intense precursor). Value 1 uses
#' 			the 2nd-most intense precursor for a chosen ion (and its data-dependent scans)
#' 			, etc.
#' @param useRtLimit Whether to enforce the given retention time window.
#' @param Args A list of arguments that will be handed to the xcms-method findPeaks via do.call
#' @param settings Options to be used for processing. Defaults to the options loaded via
#' 			\code{\link{loadRmbSettings}} et al. Refer to there for specific settings.
#' @param progressbar The progress bar callback to use. Only needed for specialized applications.
#' 			Cf. the documentation of \code{\link{progressBarHook}} for usage.
#' @return The \code{msmsWorkspace} with msms-spectra read.
#' @seealso \code{\link{msmsWorkspace-class}}, \code{\link{msmsWorkflow}}
#' @author Michael Stravs, Eawag <michael.stravs@@eawag.ch>
#' @author Erik Mueller, UFZ
#' @export
msmsRead <- function(w, filetable = NULL, files = NULL, cpdids = NULL, 
					readMethod, mode, confirmMode = FALSE, useRtLimit = TRUE, 
					Args = NULL, settings = getOption("RMassBank"), progressbar = "progressBarHook"){
	
	##Read the files and cpdids according to the definition
	##All cases are silently accepted, as long as they can be handled according to one definition
	if(is.null(filetable)){
		##If no filetable is supplied, filenames must be named explicitly
		if(is.null(files))
			stop("Please supply the files")
		
		##Assign the filenames to the workspace
		w@files <- unlist(files)
		
		##If no filetable is supplied, cpdids must be delivered explicitly or implicitly within the filenames
		if(is.null(cpdids)){
			splitfn <- strsplit(files,"_")
			splitsfn <- sapply(splitfn, function(x) x[length(x)-1])
			if(suppressWarnings(!is.na(as.numeric(splitsfn))))
				stop("Please supply the cpdids corresponding to the files in the filetable or the filenames")
			cpdids <- splitsfn
		}
	} else{
		##If a filetable is supplied read it
		tab <- read.csv(filetable)
		w@files <- tab[,"Files"]
		cpdids <- tab[,"ID"]
	}
	
	##If there's more cpdids than filenames or the other way around, then abort
	if(length(w@files) != length(cpdids)){
		stop("There are a different number of cpdids than files")
	}
	if(!(readMethod %in% c("mzR","peaklist","xcms"))){
		stop("The supplied method does not exist")
	}
	if(!all(file.exists(w@files))){
		stop("The supplied files don't exist")
	}
	
	##Progressbar
	nLen <- length(w@files)
	nProg <- 0
	pb <- do.call(progressbar, list(object=NULL, value=0, min=0, max=nLen))

	##This should work
	if(readMethod == "mzR"){
		count <- 1
		w@specs <-  lapply(w@files, function(fileName){
			spec <- findMsMsHR(fileName, cpdids[count], mode, confirmMode, useRtLimit,
		 		ppmFine = settings$findMsMsRawSettings$ppmFine,
		 		mzCoarse = settings$findMsMsRawSettings$mzCoarse,
		 		fillPrecursorScan = settings$findMsMsRawSettings$fillPrecursorScan,
		 		rtMargin = settings$rtMargin,
		 		deprofile = settings$deprofile)
			
			## Progress:
			nProg <<- nProg + 1
			pb <- do.call(progressbar, list(object=pb, value= nProg))
			
			##Counting the index of cpdids
			count <<- count + 1
			return(spec)
		})
		names(w@specs) <- basename(as.character(w@files))
		return(w)
	}
	
	##Peaklist-readmethod 
	if(readMethod == "peaklist"){
		count <- 1
		w@specs <-  lapply(w@files, function(fileName){
			spec <- addpeaksmanually(files, cpdids[count], mode=mode, findPeaksArgs=Args)
			count <<- count + 1
			return(spec)
		})
		return(w)
	}
	
	##xcms-readmethod 
	if(readMethod == "xcms"){
		count <- 1
		w@specs <-  lapply(w@files, function(fileName){
			spec <- findMsMsHRperxcms.direct(fileName, cpdids[count], mode=mode, findPeaksArgs=Args)
			count <<- count + 1
			return(spec)
		})
		return(w)
	}
}