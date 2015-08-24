# TODO: Add comment
# 
# Author: stravsmi
###############################################################################


#' Combine workspaces for multiplicity filtering
#' 
#' Combines multiple msmsWorkspace items to one workspace which is used for
#' multiplicity filtering. 
#' 
#' This feature is particularily meant to be used in 
#' conjunction with the \code{confirmMode} option of \code{\link{msmsWorkflow}}:
#' a file can be analyzed with \code{confirmMode = 0} (default) and subsequently
#' with \code{confirmMode = 1} (take second highest scan). The second analysis
#' should contain "the same" spectra as the first one (but less intense) and can
#' be used to confirm the peaks in the first spectra.  
#' 
#' TO DO: Enable the combination of workspaces for combining e.g. multiple
#' energy settings measured separately.
#' 
#' @usage combineMultiplicities(workspaces)  
#' @param workspaces A vector of \code{msmsWorkspace} items. The first item is
#' 		taken as the "authoritative" workspace, i.e. the one which will be used
#' 		for the record generation. The subsequent workspaces will only be used
#' 		for multiplicity filtering.
#' @seealso \code{\link{msmsWorkspace-class}}
#' @return A \code{msmsWorkspace} object prepared for step 8 processing.
#' 
#' @examples \dontrun{
#' 	w <- newMsmsWorkspace
#'  w@@files <- c("spec1", "spec2")
#'  w1 <- msmsWorkflow(w, steps=c(1:7), mode="pH")
#'  w2 <- msmsWorkflow(w, steps=c(1:7), mode="pH", confirmMode = 1)
#'  wTotal <- combineMultiplicities(c(w1, w2))
#'  wTotal <- msmsWorkflow(wTotal, steps=8, mode="pH", archivename = "output")
#'  # continue here with mbWorkflow 
#' }
#' 
#' @author Stravs MA, Eawag <michael.stravs@@eawag.ch>
#' @export
combineMultiplicities <- function(workspaces)
{
	# get the first workspace for output
	wOut <- workspaces[[1]]
	# Combine: aggregated peaks
	wOut@aggregated = do.call(rbind,
			lapply(workspaces, function(w) w@aggregated))
	
	return(wOut)
}


#' Determine processed steps
#' 
#' This function reads out the content of different slots of the \code{workspace}
#' object and finds out which steps have already been processed on it.
#' 
#' @param workspace A \code{msmsWorkspace} object. 
#' 
#' @return An array containing all \code{msmsWorkflow} steps which have 
#' likely been processed.  
#' 
#' @examples \dontrun{
#' findProgress(w)
#' }
#' @author Stravs MA, Eawag <michael.stravs@@eawag.ch>
#' @export
findProgress <- function(workspace)
{
	step4 <- (!is.null(workspace@parent))
	
	if(step4)
		w123 <- workspace@parent
	else
		w123 <- workspace
    step1 <- (length(w123@spectra) > 0)
    step2 <- step1 && !all(is.na(lapply(w123@spectra, function(s) s@empty))) 
    step3 <- (nrow(w123@aggregated) > 0)
	
	step5 <- step4 && !all(is.na(lapply(workspace@spectra, function(s) s@empty))) 
	step6 <- step5 && (nrow(workspace@aggregated) > 0)
	step7 <- step6 && ("matchedReanalysis" %in% colnames(workspace@aggregated))
	step8 <- step7 && ("filterOK" %in% colnames(workspace@aggregated))
	
    steps <- which(c(step1, step2, step3, step4, step5, step6, step7, step8))
    return(steps)
}

#' Update settings to current version
#'
#' Checks if all necessary fields are present in the current settings
#' and fills in default values from the \code{\link{RmbDefaultSettings}}
#' if required.
#' 
#' @note Important: There is a change in behaviour of RMassBank in certain cases when \code{filterSettings} is not
#' present in the old settings! The default pre-recalibration cutoff from \code{\link{RmbDefaultSettings}} is 10000.
#' Formerly the pre-recalibration cutoff was set to be 10000 for positive spectra but 0 for negative spectra.
#' 
#' Updating the settings files is preferred to using the \code{updateSettings} function.
#' 
#' @param settings The set of settings to check and update.
#' 
#' @param warn Whether to update parameters quietly (\code{FALSE}) or to notify the user
#' 	of the changed parameters (\code{TRUE}, default.) This serves to make the user aware that
#' standard parameters are filled in!
#' 
#' @return The updated set of settings.
#' 
#' @examples \dontrun{
#' w@@settings <- updateSettings(w@@settings)
#' }
#' 
#' @author Stravs MA, Eawag <michael.stravs@@eawag.ch>
#' @export
#' 
updateSettings <- function(settings, warn=TRUE)
{
  settings.new <- .settingsList
  settings.old <- settings
  renew <- setdiff(names(settings.new), names(settings.old))
  if(length(renew) > 0 && warn==TRUE){
    warning(paste0("Your settings are outdated! The following fields were taken from default values: ", 
            paste(renew,collapse=", ")))
    if("filterSettings" %in% renew)
      warning("The default values of filterSettings could change the processing behaviour if you have negative-mode spectra. Check ?updateSettings for details.")
  }
  settings.old[renew] <- settings.new[renew]
  return(settings.old)
}
