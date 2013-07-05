

#' Workspace for \code{msmsWorkflow} data
#' 
#' A workspace which stores input and output data for \code{\link{msmsWorkflow}}.
#' 
#' Slots:
#' 
#'  \describe{
#' 	\item{files}{The input file names}
#' 	\item{specs}{The spectra extracted from the raw files}
#'  \item{analyzedSpecs}{The spectra with annotated peaks after workflow step 2.}
#'  \item{aggregatedSpecs}{The \code{analyzedSpec} data regrouped and aggregated, 
#' 		after workflow step 3.}
#'  \item{rc, rc.ms1}{The recalibration curves generated in workflow step 4.}
#'  \item{recalibratedSpecs}{The spectra from \code{specs} recalibrated with the curves
#' 		from \code{rc, rc,ms1}.} 
#' 	\item{analyzedRcSpecs}{The recalibrated spectra with annotated peaks after 
#' 		workflow step 5.}
#' 	\item{aggregatedRcSpecs}{The \code{analyzedRcSpec} data regrouped and aggregated, 
#' 		after workflow step 6.}
#'  \item{reanalyzedRcSpecs}{The regrouped and aggregated spectra, with added reanalyzed
#' 		peaks (after step 7, see \code{\link{reanalyzeFailpeaks}}).}
#'  \item{refilteredRcSpecs}{Final data to use for MassBank record creation after 
#' 		multiplicity filtering (step 8).}
#' }
#' 
#' Methods: \describe{
#' 	\item{show}{Shows a brief summary of the object. Currently only the included files.}
#' 	}
#' 
## ' @slot files The input file names
## ' @slot specs The spectra extracted from the raw files
## ' @slot analyzedSpecs The spectra with annotated peaks after workflow step 2.
## ' @slot aggregatedSpecs The \code{analyzedSpec} data regrouped and aggregated, 
## ' 		after workflow step 3.
## ' @slot rc,rc.ms1 The recalibration curves generated in workflow step 4.
## ' @slot recalibratedSpecs The spectra from \code{specs} recalibrated with the curves
## ' 		from \code{rc, rc,ms1}.
## ' @slot analyzedRcSpecs The recalibrated spectra with annotated peaks after 
## ' 		workflow step 5.
## ' @slot aggregatedRcSpecs The \code{analyzedRcSpec} data regrouped and aggregated, 
## ' 		after workflow step 6.
## ' @slot reanalyzedRcSpecs The regrouped and aggregated spectra, with added reanalyzed
## ' 		peaks (after step 7, see \code{\link{reanalyzeFailpeaks}}).
## ' @slot refilteredRcSpecs Final data to use for MassBank record creation after 
## ' 		multiplicity filtering (step 8).
#' 
## ' @method show,msmsWorkflow Shows a brief summary of the object. Currently only the included files.
#' 
#' @seealso \code{\link{msmsWorkflow}}
#' @author Michael Stravs, Eawag <michael.stravs@@eawag.ch>
#' @name msmsWorkspace-class
#' @docType class
#' @exportClass msmsWorkspace
#' @export
setClass("msmsWorkspace",
		representation(
				files = "character",
				specs = "list",
				analyzedSpecs = "list",
				aggregatedSpecs = "list",
				rc = "ANY",
				rc.ms1 = "ANY",
				recalibratedSpecs = "list",
				analyzedRcSpecs = "list",
				aggregatedRcSpecs = "list",
				reanalyzedRcSpecs = "list",
				refilteredRcSpecs = "list",
				archivename = "character",
				settings = "list"
				),
		)

#' Workspace for \code{mbWorkflow} data
#' 
#' A workspace which stores input and output data for use with \code{mbWorkflow}.
#' 
#' Slots:
#'  \describe{
#' 	\item{aggregatedRcSpecs, refilteredRcSpecs}{The corresponding
#' 		 input data from \code{\link{msmsWorkspace-class}}}
#'  \item{additionalPeaks}{A list of additional peaks which can be loaded
#' 		using \code{\link{addPeaks}}.}
#'  \item{mbdata, mbdata_archive, mbdata_relisted}{Infolist data: Data for
#' 		annotation of MassBank records, which can be loaded using
#' 		\code{\link{loadInfolists}}.}
#'  \item{compiled, compiled_ok}{
#' 		Compiled tree-structured MassBank records. \code{compiled_ok} contains
#' 		only the compounds with at least one valid spectrum.}
#'  \item{mbfiles}{Compiled MassBank records in text representation.}
#'  \item{molfile}{MOL files with the compound structures.}
#'  \item{ok,problems}{Index lists for internal use which denote which compounds
#' 		have valid spectra.}
#' }
#' 
## ' @method show,mbWorkflow Shows a brief summary of the object. Currently only a stub. 
#' Methods:
#' \describe{
#' 		\item{show}{Shows a brief summary of the object. Currently only a stub.}
#' }
#' 
#' @seealso \code{\link{mbWorkflow}}
#' 
#' @docType class
#' @author Michael Stravs, Eawag <michael.stravs@@eawag.ch>
#' @name mbWorkspace-class
#' @exportClass mbWorkspace
setClass("mbWorkspace",
		representation(
				# input data:
				aggregatedRcSpecs = "list",
				refilteredRcSpecs = "list",
				additionalPeaks = "data.frame", # ex additional_peaks
				# infolists data:
				mbdata = "list",
				mbdata_archive = "data.frame",
				mbdata_relisted = "list",				
				# output data:
				compiled = "list",
				compiled_ok = "list",
				mbfiles = "list",
				molfile = "list",
				ok = "integer",
				problems = "integer"
				)
		)
		

#' Create new empty workspace or load saved data for \code{msmsWorkflow}
#' 
#' Creates an empty workspace or loads an existing workspace from disk.
#' 
#' \code{newMsmsWorkspace} creates a new empty workspace for use with \code{msmsWorkflow.}
#' 
#' \code{loadMsmsWorkspace} loads a workspace saved using \code{\link{archiveResults}}. 
#' Note that it also successfully loads data saved with the old RMassBank data format
#' into the new \code{msmsWorkspace} object.   
#' 
#' @aliases newMsmsWorkspace loadMsmsWorkspace
#' @param files If given, the files list to initialize the workspace with.
#' @return A new \code{msmsWorkspace} object
#' @seealso \code{\link{msmsWorkflow}}, \code{\link{msmsWorkspace-class}}
#' 
#' @author Michael Stravs, Eawag <michael.stravs@@eawag.ch>
#' @export
newMsmsWorkspace <- function(files = character(0))
{
	w <- new("msmsWorkspace")
	w@files <- files
	return(w)
}

#' @export
loadMsmsWorkspace <- function(fileName, loadSettings = FALSE)
{
	tempEnv <- new.env()
	load(fileName, envir=tempEnv)
	# Look if there is a msmsWorkspace in the file
	objs <- ls(tempEnv)
	isWs <- unlist(lapply(objs, function(obj) "msmsWorkspace" %in% class(tempEnv[[obj]])))
	whichWs <- match(TRUE, isWs)
	# Found? Then just return it.
	if(!is.na(whichWs))
	{
		w <- tempEnv[[objs[[whichWs]]]]
		# If settings don't have to be loaded, erase them from the object
		if(loadSettings == FALSE)
			w@settings <- list()
	}
	# Otherwise hope to load the dataset into a new workspace
	else
	{
		w <- new("msmsWorkspace")
		dataset <- c("rc",
				"rc.ms1",
				"files",
				"specs",
				"analyzedSpecs",
				"aggregatedSpecs",
				"recalibratedSpecs",
				"analyzedRcSpecs",
				"aggregatedRcSpecs",
				"reanalyzedRcSpecs",
				"refilteredRcSpecs",
				"archivename")
		for(var in dataset)
		{
			if(exists(var, envir=tempEnv))
				slot(w, var) <- tempEnv[[var]]
		}
		# Check if settings exist...
		if((loadSettings == TRUE) && exists("RmbSettings", envir=tempEnv))
			w@settings <- tempEnv$RmbSettings
	}
	# If loadSettings is set: load the settings into RMassBank
	if((loadSettings == TRUE) && (length(w@settings) > 0))
		loadRmbSettings(w@settings)
	else if (loadSettings == TRUE)
		warning("You specified to load RMassBank settings from the saved workspace, but no stored settings were found.")
	return(w)
}

#' Create new workspace for \code{mbWorkflow}
#' 
#' Creates a new workspace for use with \code{\link{mbWorkflow}}.
#' 
#' The workspace input data will be  loaded from the \code{\link{msmsWorkspace-class}}
#' object provided by the parameter \code{w}.
#' 
#' @param w The input \code{msmsWorkspace} to load input data from.
#' @return A new \code{mbWorkflow} object with the loaded input data.
#' 
#' @seealso \code{\link{mbWorkflow}}, \code{\link{msmsWorkspace-class}}
#' @author Michael Stravs, Eawag <michael.stravs@@eawag.ch>
#' @export
newMbWorkspace <- function(w)
{
	mb <- new("mbWorkspace",
			aggregatedRcSpecs = w@aggregatedRcSpecs,
			refilteredRcSpecs = w@refilteredRcSpecs)
	return(mb)
}


#' @name show,msmsWorkspace-method
#' @aliases show,msmsWorkspace-method
#' @docType methods
#' @rdname msmsWorkspace-class
#' @export
setMethod("show", "msmsWorkspace",
		function(object) {
			cat("Object of class \"",class(object),"\"\n",sep="")
			cat(" with files: \n")
			sapply(basename(object@files), function(x) cat("  -", x, "\n"))
		})  



#'
#' @name show,mbWorkspace-method 
#' @aliases show,mbWorkspace-method
#' @rdname mbWorkspace-class
#' @docType methods
#' @export
setMethod("show", "mbWorkspace",
		function(object) {
			cat("Object of class \"",class(object),"\"\n",sep="")
		}) 
		
