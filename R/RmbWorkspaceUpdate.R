# TODO: Add comment
# 
# Author: stravsmi
###############################################################################


.updateObject.RmbWorkspace <- setMethod("updateObject", signature(object="msmsWorkspace"), function(object, ..., verbose = FALSE) 
		{
			w <- object
			if(isVersioned(w))
				if(all(isCurrent(w)))
					return(w)
			# get msmsWorkspace version
			if(!isVersioned(w))
				v <- "1.0.0"
			else
				v <- classVersion(w)["msmsWorkspace"]
			w.new <- w
			# gradually step up versions
			# 2.0.1: 
			# * spectra go from specs, analyzedSpecs or their rc analogs to Spectra
			# * data pre recalibration get shifted to "parent workspace"
			if(v < "2.0.1")
			{
				w.old <- w.new
				w.new <- new("msmsWorkspace")
				w.new@files <- w.old@files
				# Do we have recalibration done? If so: all data in the WS will be the recalibrated data, the unrecalibrated data will be
				# moved into a new parent workspace which is referenced
				progress <- .findProgress.v1(w.old)
				if(4 %in% progress)
				{
					w.parent <- w.old
					w.parent@recalibratedSpecs <- list()
					w.parent@analyzedRcSpecs <- list()
					w.parent@aggregatedRcSpecs <- list()
					w.parent@reanalyzedRcSpecs <- list()
					w.parent@refilteredRcSpecs <- list()
					w.old@specs <- w.old@recalibratedSpecs
					w.old@analyzedSpecs <- w.old@analyzedRcSpecs
					w.parent.new <- updateObject(w.parent)
					w.new@parent <- w.parent.new
				}
				w.new@spectra <- .updateObject.spectra(w.old@specs, w.old@analyzedSpecs)
			}
			
			return(w.new)
		})



.updateObject.spectra <- function(specs, analyzedSpecs)
{
	if((length(specs) != length(analyzedSpecs)) && (0 != length(analyzedSpecs) ))
		stop("updateObject: Could not update object because data is inconsistent. length(analyzedSpecs) != length(specs) or 0")
	# process info ex specs
	spectra <- lapply(specs, function(spec){
				set <- new("RmbSpectraSet")
				# identifiers and properties
				set@mz <- spec$mz$mzCenter
				set@id <- as.integer(spec$id)
				set@formula <- spec$formula
				set@found <- as.logical(spec$foundOK)
				# now parent and child MS
				# check for parent recalibration column
				if("mzRecal" %in% colnames(spec$parentPeak))
					mzcol <- "mzRecal"
				else
					mzcol <- "mz"
				set@parent <- new("Spectrum1", 
						mz = spec$parentPeak[,mzcol],
						intensity = spec$parentPeak[,2],
						polarity = as.integer(spec$parentHeader$polarity),
						peaksCount = as.integer(spec$parentHeader$peaksCount),
						rt = spec$parentHeader$retentionTime,
						acquisitionNum = as.integer(spec$parentHeader$acquisitionNum),
						tic = spec$parentHeader$totIonCurrent,
						centroided = TRUE
						)
				# get MSMS data from spec$peaks into RmbSpectrum2 objects
				children.p1 <- lapply(spec$peaks, function(peaks)
						{
							if("mzRecal" %in% colnames(peaks))
								mzcol <- "mzRecal"
							else
								mzcol <- "mz"
							new("RmbSpectrum2",
									mz=peaks[,mzcol],
									intensity=peaks[,2],
									peaksCount=nrow(peaks))
						})
				# get header data from spec$childHeaders into separate RmbSpectrum2 objects
				children.p2 <- apply(spec$childHeaders, 1, function(line)
						{
							new("RmbSpectrum2",
									precScanNum = as.integer(line["precursorScanNum"]),
									precursorMz = line["precursorMZ"],
									precursorIntensity = line["precursorIntensity"],
									precursorCharge = as.integer(line["precursorCharge"]),
									collisionEnergy = line["collisionEnergy"],
									tic = line["totIonCurrent"],
									rt = line["retentionTime"],
									acquisitionNum = as.integer(line["acquisitionNum"]),
									centroided = TRUE
							)
						}) 
				# merge MSMS RmbSpectrum2 with header RmbSpectrum2
				children <- mapply(function(c1,c2)
						{
							c2slots <- c("precScanNum","precursorMz", "precursorIntensity", "precursorCharge", "collisionEnergy",
									"tic", "rt", "acquisitionNum", "centroided")
							for(c2slot in c2slots)
								slot(c1, c2slot) <- slot(c2, c2slot)
							return(c1)
						}, children.p1, children.p2)
				set@children <- as(children, "SimpleList")
				
				return(set)
			})
	spectra <- mapply(function(set, name)
			{
				set@name <- name
				return(set)
			},
			spectra, names(specs))
	
				
	# add info ex analyzedSpecs if present
	if(length(analyzedSpecs) > 0)
	{
		spectra <- mapply(function(set, analyzedSpec)
				{
					if(length(analyzedSpec$msmsdata) != length(set@children))
						stop("updateObject: Could not update object because data is inconsistent. length(analyzedSpec$msmsdata) != length(set@children)")
					children <- mapply(function(spectrum, msmsrecord)
							{
								if(!is.data.frame(msmsrecord$childBad))
									msmsrecord$childBad <- data.frame()
								# note: mz/intensity are replaced with the values from the analyzed spectrum,
								# such as to have a mass multiple times for multiple matched formulas
								
								# check if the spectrum has recalibrated masses; if yes, use those
								if("mzRecal" %in% colnames(msmsrecord$childRaw))
									mzcol <- "mzRecal"
								else
									mzcol <- "mz"
								spectrum@mz <- c(msmsrecord$childFilt$mzFound,
										msmsrecord$childBad$mzFound,
										msmsrecord$childUnmatched$mzFound,
										msmsrecord$childRawLow[,mzcol], 
										msmsrecord$childRawSatellite[,mzcol])
								spectrum@intensity <- c(msmsrecord$childFilt$int,
										msmsrecord$childBad$int,
										msmsrecord$childUnmatched$int,
										msmsrecord$childRawLow$int, 
										msmsrecord$childRawSatellite$int)
								spectrum@peaksCount <- length(spectrum@mz)

								spectrum@satellite <- as.logical(c(
												rep(FALSE,nrow(msmsrecord$childFilt)),
												rep(FALSE,nrow(msmsrecord$childBad)),
												rep(FALSE,nrow(msmsrecord$childUnmatched)),
												rep(FALSE,nrow(msmsrecord$childRawLow)), 
												rep(TRUE,nrow(msmsrecord$childRawSatellite)))
								)
								
								
								spectrum@low <- as.logical(c(
												rep(FALSE,nrow(msmsrecord$childFilt)),
												rep(FALSE,nrow(msmsrecord$childBad)),
												rep(FALSE,nrow(msmsrecord$childUnmatched)),
												rep(TRUE,nrow(msmsrecord$childRawLow)), 
												rep(FALSE,nrow(msmsrecord$childRawSatellite)))
								)
								
								spectrum@rawOK <- as.logical(c(
												rep(TRUE,nrow(msmsrecord$childFilt)),
												rep(TRUE,nrow(msmsrecord$childBad)),
												rep(TRUE,nrow(msmsrecord$childUnmatched)),
												rep(FALSE,nrow(msmsrecord$childRawLow)), 
												rep(FALSE,nrow(msmsrecord$childRawSatellite)))
								)
								
								spectrum@good <- as.logical(c(
										msmsrecord$childFilt$good,
										msmsrecord$childBad$good,
										msmsrecord$childUnmatched$good,
										rep(NA,nrow(msmsrecord$childRawLow)), 
										rep(NA,nrow(msmsrecord$childRawSatellite)))
								)
								spectrum@mzCalc <- as.numeric(c(
										msmsrecord$childFilt$mzCalc,
										msmsrecord$childBad$mzCalc,
										msmsrecord$childUnmatched$mzCalc,
										rep(NA,nrow(msmsrecord$childRawLow)), 
										rep(NA,nrow(msmsrecord$childRawSatellite)))
								)
								spectrum@formula <- as.character(c(
										msmsrecord$childFilt$formula,
										msmsrecord$childBad$formula,
										msmsrecord$childUnmatched$formula,
										rep(NA,nrow(msmsrecord$childRawLow)), 
										rep(NA,nrow(msmsrecord$childRawSatellite)))
								)
								
								spectrum@dbe <- as.numeric(c(
												msmsrecord$childFilt$dbe,
												msmsrecord$childBad$dbe,
												msmsrecord$childUnmatched$dbe,
												rep(NA,nrow(msmsrecord$childRawLow)), 
												rep(NA,nrow(msmsrecord$childRawSatellite)))
								)
								
								spectrum@formulaCount <- as.integer(c(
												msmsrecord$childFilt$formulaCount,
												msmsrecord$childBad$formulaCount,
												msmsrecord$childUnmatched$formulaCount,
												rep(NA,nrow(msmsrecord$childRawLow)), 
												rep(NA,nrow(msmsrecord$childRawSatellite)))
								)
								
								spectrum@dppm <- as.numeric(c(
												msmsrecord$childFilt$dppm,
												msmsrecord$childBad$dppm,
												msmsrecord$childUnmatched$dppm,
												rep(NA,nrow(msmsrecord$childRawLow)), 
												rep(NA,nrow(msmsrecord$childRawSatellite)))
								)
								spectrum@dppmBest <- as.numeric(c(
												msmsrecord$childFilt$dppmBest,
												msmsrecord$childBad$dppmBest,
												msmsrecord$childUnmatched$dppmBest,
												rep(NA,nrow(msmsrecord$childRawLow)), 
												rep(NA,nrow(msmsrecord$childRawSatellite)))
								)
								
#								.RmbSpectrum2 <- setClass("RmbSpectrum2",
#										representation = representation(
								## satellite="logical",
								## low="logical",
								## rawOK ="logical",
								## good = "logical",
								## mzCalc = "numeric",
								## formula = "character",
								## formulaCount = "integer",
								## dppm = "numeric",
								## dppmBest = "numeric",
#										),
								
								return(spectrum)
								
							},
							set@children, analyzedSpec$msmsdata)
							set@children <- as(children, "SimpleList")
							return(set)
				},
				spectra, analyzedSpecs)
	}


	return(as(spectra, "SimpleList"))
	
}

# Finds progress in the "old workspace version" to determine whether to take the old spectra or the recalibrated ones (and
# make a parent workspace)
.findProgress.v1 <- function(workspace)
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