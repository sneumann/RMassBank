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
#'  wTotal <- msmsWorkflow(wTotal, steps=8, mode="pH", archiveName = "output")
#'  # continue here with mbWorkflow 
#' }
#' 
#' @author Stravs MA, Eawag <michael.stravs@@eawag.ch>
#' @export
combineMultiplicities <- function(workspaces)
{
	# get the first workspace for output
	wOut <- workspaces[[1]]
	# Combine: reanalyzedRcSpecs$peaksMatched, peaksMatchedReanalysis, peaksReanalyzed
	# (which are the result of step 7 and used in step 8)
	wOut@reanalyzedRcSpecs$peaksMatched = do.call(rbind,
			lapply(workspaces, function(w) w@reanalyzedRcSpecs$peaksMatched))
	wOut@reanalyzedRcSpecs$peaksMatchedReanalysis = do.call(rbind,
			lapply(workspaces, function(w) w@reanalyzedRcSpecs$peaksMatchedReanalysis))
	wOut@reanalyzedRcSpecs$peaksReanalyzed = do.call(rbind,
			lapply(workspaces, function(w) w@reanalyzedRcSpecs$peaksReanalyzed))
	
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
#' @author Stravs MA, Eawag <michael.stravs@@eawag.ch>
#' @export
findProgress <- function(workspace)
{
    step1 <- (length(workspace@specs) > 0)
    step2 <- (length(workspace@analyzedSpecs) > 0)
    step3 <- (length(workspace@aggregatedSpecs) > 0)
    step4 <- (length(workspace@recalibratedSpecs) > 0)
    step5 <- (length(workspace@analyzedRcSpecs) > 0)
    step6 <- (length(workspace@aggregatedRcSpecs) > 0)
    step7 <- (length(workspace@reanalyzedRcSpecs) > 0)
    step8 <- (length(workspace@refilteredRcSpecs) > 0)
    steps <- which(c(step1, step2, step3, step4, step5, step6, step7, step8))
    return(steps)
}
