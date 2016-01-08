# Based on the code from MSnbase
# * contributed by Guangchuang Yu <guangchuangyu@gmail.com>
# * Modified by Sebastian Gibb <mail@sebastiangibb.de>
# * adapted into RMassBank by Michele Stravs <stravsmi@eawag.ch>
## setMethod("writeMgfData",
##           signature = signature("Spectrum1"),
##           function(object,
##                    con = "spectrum.mgf",
##                    COM = NULL,
##                    TITLE = NULL) {
##             writeMgfDataFile1(list(object), con = con, COM = COM, TITLE = TITLE,
##                              verbose = FALSE)
##           })
## 
  
setMethod("writeMgfData",
		signature=signature("RmbSpectraSet"),
		function(object,
				con="spectrum.mgf",
				COM = NULL,
				TITLE = NULL) {
			writeMgfSpectraSet(object,
					con = con, COM = COM, TITLE = TITLE,
					verbose = FALSE)
		})

writeMgfSpectraSet <- function(object,
					con = con, COM = COM, TITLE = TITLE,
					verbose = FALSE)
{
	if (class(con) == "character" && file.exists(con)) {
		message("Overwriting ", con, "!")
		unlink(con)
	}
	
	con <- file(description = con, open = "at")
	on.exit(close(con))
	
	if (is.null(COM)) {
		COM <- paste0("COM=", object@id, " ", object@name, " ", object@mz, " ", object@formula,
				"exported by MSnbase/RMassBank on ", date())
	}
	cat(COM, file = con, sep = "")
	writeMgfContent1(object@parent, TITLE = NULL, con = con)
	for(chi in as.list(object@children))
	{
		MSnbase:::writeMgfContent(chi, TITLE=NULL, con=con)
	}
	
	
}
## 
## writeMgfDataFile1 <- function(splist, con, COM = NULL, TITLE = NULL,
##                              verbose = TRUE) {
##   if (class(con) == "character" && file.exists(con)) {
##     message("Overwriting ", con, "!")
##     unlink(con)
##   }
## 
##   con <- file(description = con, open = "at")
##   on.exit(close(con))
## 
##   if (is.null(COM)) {
##     COM <- paste0("COM=", ifelse(length(splist) <= 1, "Spectrum", "Experiment"),
##                   "exported by MSnbase on ", date())
##   }
##   cat(COM, file = con, sep = "")
## 
##   verbose <- verbose & length(splist) > 1
## 
##   if (verbose)
##     pb <- txtProgressBar(min = 0, max = length(splist), style = 3)
## 
##   for (i in seq(along=splist)) {
##     if (verbose)
##       setTxtProgressBar(pb, i)
## 
##     writeMgfContent1(splist[[i]], TITLE = NULL, con = con)
##   }
##   if (verbose)
##     close(pb)
## }
## 
## writeMgfContent1 <- function(sp, TITLE = NULL, con) {
##   .cat <- function(..., file=con, sep="", append=TRUE) {
##     cat(..., file=file, sep=sep, append=append)
##   }
## 
##   .cat("\nBEGIN IONS\n",
##        "SCANS=", acquisitionNum(sp))
## 
##   if (is.null(TITLE)) {
##     .cat("\nTITLE=msLevel ", msLevel(sp),
##          "; retentionTime ", rtime(sp),
##          "; scanNum ", acquisitionNum(sp))
## 
##     if (length(scanIndex(sp))) {
##       .cat("; scanIndex ", scanIndex(sp))
##     }
## 
##     if (msLevel(sp) > 1) {
##       .cat("; precMz ", precursorMz(sp),
##            "; precCharge ", precursorCharge(sp))
##     }
##   } else {
##     .cat("\nTITLE=", TITLE)
##   }
## 
##   .cat("\nRTINSECONDS=", rtime(sp))
## 
## 
##   .cat("\n", paste(mz(sp), intensity(sp), collapse = "\n"))
##   .cat("\nEND IONS\n")
## }
