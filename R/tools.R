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
