# TODO: Add comment
# 
# Author: stravsmi
###############################################################################


#' Select peaks from aggregate table
#' 
#' Selects peaks from aggregate table according to different criteria.
#' 
#' @param o \code{msmsWorkspace} or aggregate \code{data.frame} from a workspace.
#' @param good if \code{TRUE}, include good (matched within filter criteria) peaks.
#' @param bad if \code{TRUE}, include bad (not matched within filter criteria) peaks. Note: \code{good} and \code{bad} can be
#' 			combined, both are returned in that case.
#' @param cleaned if \code{TRUE}, return only peaks which passed the noise filter. Note: If the noise filter was not applied, the
#' 			parameter has no effect. Also, a \code{noise} column is in any case added to the output, even if not present before.
#' @param best if \code{TRUE}, only select the best match for each peak (i.e. the formula with smallest delta ppm). Otherwise multiple
#' 			matches can be returned.
#' @param ... no additional parameters
#' @return Peak dataframe according to the specified criteria.
#' 
#' @author stravsmi
#' @export
setGeneric("selectPeaks", function(o, ...) standardGeneric("selectPeaks"))

#' Add and initialize dataframe column
#' 
#' Adds a new column of a defined type to a \code{data.frame} and initializes it to a value.
#' The advantage of doing this over adding it with \code{$} or \code{[,""]} is that the case
#' \code{nrow(o) == 0} is adequately handled and doesn't raise an error. 
#' 
#' @param o \code{data.frame} to add the column to
#' @param name Name of the new column
#' @param type Data type of the new column
#' @param value Initial value of the new column (\code{NA} if not given)
#' @return Expanded data frame.
#' 
#' @author stravsmi
#' @export
setGeneric("addProperty", function(o, name, type, value=NA) standardGeneric("addProperty"))

setGeneric("property", function(o, property) standardGeneric("property"))
setGeneric("property<-", function(o, property, value, addNew = FALSE, class="") standardGeneric("property<-"))