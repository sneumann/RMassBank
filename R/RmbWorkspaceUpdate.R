# TODO: Add comment
# 
# Author: stravsmi
###############################################################################


.updateObject.RmbWorkspace <- setMethod("updateObject", signature(object="msmsWorkspace"), function(object, ..., verbose = FALSE) 
		{
			w <- object
			w <- .updateObject.RmbWorkspace.1to2(w, ..., verbose)
			w <- .updateObject.RmbWorkspace.RmbSpectraSet(w, ..., verbose)
			classVersion(w)["msmsWorkspace"] <- "2.0.4"
			w
		})


.updateObject.RmbWorkspace.RmbSpectraSet <- function(object, ..., verbose = FALSE) 
{
	w <- object
	v <- classVersion(w)["msmsWorkspace"]
	if(v < "2.0.4")
	{
		for(i in seq_len(length(w@spectra)))
			w@spectra[[i]] <- updateObject(w@spectra[[i]])
	}
	w
}

.updateObject.RmbWorkspace.1to2 <- function(object, ..., verbose = FALSE) 
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
			# remove the recalibrated processing results from the future parent workspace
			slot(w.parent, "recalibratedSpecs", check=FALSE) <- NULL
			slot(w.parent, "analyzedRcSpecs", check=FALSE) <- NULL
			slot(w.parent, "aggregatedRcSpecs", check=FALSE) <- NULL
			slot(w.parent, "reanalyzedRcSpecs", check=FALSE) <- NULL
			slot(w.parent, "refilteredRcSpecs", check=FALSE) <- NULL
			# move the recalibrated base data to the base data of the current workspace
			slot(w.old, "specs", check=FALSE) <- w.old@recalibratedSpecs
			slot(w.old, "analyzedSpecs", check=FALSE) <- w.old@analyzedRcSpecs
			w.parent.new <- updateObject(w.parent)
			# fill in the calibrations into the parent, which are otherwise not copied
			slot(w.parent.new, "rc", check=FALSE) <- w.old@rc
			slot(w.parent.new, "rc.ms1", check=FALSE) <- w.old@rc.ms1
			w.new@parent <- w.parent.new
		}
		w.new@spectra <- .updateObject.spectra(w.old@specs, w.old@analyzedSpecs)
		if(any(c(3,6,7) %in% progress))
		{
			w.new@aggregated <- aggregateSpectra(w.new, addIncomplete=TRUE)
			warning("You are loading an archive from an old RMassBank version. The aggregate tables are not loaded from the original object, but recomputed.")
			warning("If you hand-edited any aggregate table, the information might not be retained in the new object.")
		}
		# else if(6 %in% progress)
		# {
		# 	w.new@aggregated <- .updateObject.aggregated(w.old@aggregatedRcSpecs)
		# }
		# else if(3 %in% progress)
		# {
		# 	w.new@aggregated <- .updateObject.aggregated(w.old@aggregatedSpecs)
		# }
		
		if(8 %in% progress)
		{
			w.new <- .updateObject.refiltered(w.new, w.new@aggregated, w.old@refilteredRcSpecs)
			warning("You are loading an archive from an old RMassBank version. The multiplicity filtering results are not loaded from the original object, but recomputed.")
			warning("If you hand-edited any multiplicity filtering results, the information might not be retained in the new object.")
		}
		
	}
  
  
  # 2.0.4 directly from v1: update spectra polarity, because RmbSpectraSet is generated directly as a
  # 0.1.2 versioned class and does not go through the update
  w.new@spectra <- as(lapply(w.new@spectra, .updateObject.RmbSpectraSet.updatePolarity), "SimpleList")
	
	return(w.new)
}

.updateObject.spectra <- function(specs, analyzedSpecs)
{
	if((length(specs) != length(analyzedSpecs)) && (0 != length(analyzedSpecs) ))
	{
		# Try to fix this, in early processed versions it could happen that a scan got reused
		scans.analyzed <- unlist(lapply(analyzedSpecs, function(sp) sp$scan))
		scans.recorded <- unlist(lapply(specs, function(sp) sp@acquisitionNum))
		scans.reordered <- match(scans.recorded, scans.analyzed)
		analyzedSpecs <- analyzedSpecs[scans.reordered]
		id <- analyzedSpecs[[1]]$id
		warning(paste0(id, ": Spectra were reordered to match acquisition data."))
	}
	if((length(specs) != length(analyzedSpecs)) && (0 != length(analyzedSpecs) ))
	{
		stop("updateObject: Could not update object because data is inconsistent. length(analyzedSpecs) != length(specs) or 0")
		# Maybe it hasn't worked :)
	}
	
	# process info ex specs
	spectra <- lapply(specs, function(spec){
				set <- new("RmbSpectraSet")
				# identifiers and properties
				set@mz <- spec$mz$mzCenter
				set@id <- as.character(as.integer(spec$id))
				set@formula <- spec$formula
				set@found <- as.logical(spec$foundOK)
        set@smiles <- ""
				# now parent and child MS
				# check for parent recalibration column
				if(set@found)
				{
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
				}
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
					{
						# Try to fix this, in early processed versions it could happen that a scan got reused
						scans.analyzed <- unlist(lapply(analyzedSpec$msmsdata, function(sp) sp$scan))
						scans.recorded <- unlist(lapply(set@children, function(sp) sp@acquisitionNum))
						scans.reordered <- match(scans.recorded, scans.analyzed)
						analyzedSpec$msmsdata <- analyzedSpec$msmsdata[scans.reordered]
						id <- analyzedSpec$msmsdata[[1]]$id
						warning(paste0(id, ": Spectra were reordered to match acquisition data."))
					}
					
					
					if(length(analyzedSpec$msmsdata) != length(set@children))
						stop("updateObject: Could not update object because data is inconsistent. length(analyzedSpec$msmsdata) != length(set@children)")
					
					set@complete <- FALSE
					set@empty <- FALSE
					
					
					if(length(analyzedSpec$msmsdata) == 0)
					{
						empty <- TRUE
						return(set)
					}
					children <- mapply(function(spectrum, msmsrecord)
							{
								if(msmsrecord$specOK)
								
									spectrum@ok <- TRUE
								else
									spectrum@ok <- FALSE
								
								# check if the spectrum has recalibrated masses; if yes, use those
								if("mzRecal" %in% colnames(msmsrecord$childRaw))
									mzcol <- "mzRecal"
								else
									mzcol <- "mz"
								
								# create potentially missing data frames
								if(!is.data.frame(msmsrecord$childFilt))
								{
									msmsrecord$childFilt <- data.frame()
									msmsrecord$childRawLow[,"mzFound"] <- numeric()
									msmsrecord$childRawLow[,"int"] <- numeric()
								}
								if(!is.data.frame(msmsrecord$childBad))
								{
									msmsrecord$childBad <- data.frame()
									msmsrecord$childBad[,"mzFound"] <- numeric()
									msmsrecord$childBad[,"int"] <- numeric()
								}
								if(!is.data.frame(msmsrecord$childUnmatched))
								{
									msmsrecord$childUnmatched <- data.frame()
									msmsrecord$childUnmatched[,"mzFound"] <- numeric()
									msmsrecord$childUnmatched[,"int"] <- numeric()
								}
								if(!is.data.frame(msmsrecord$childRawLow))
								{
									msmsrecord$childRawLow <- data.frame()
									msmsrecord$childRawLow[,mzcol] <- numeric()
									msmsrecord$childRawLow[,"int"] <- numeric()
								}
								if(!is.data.frame(msmsrecord$childRawSatellite))
								{
									msmsrecord$childRawSatellite <- data.frame()
									msmsrecord$childRawSatellite[,mzcol] <- numeric()
									msmsrecord$childRawSatellite[,"int"] <- numeric()
									
								}
								# note: mz/intensity are replaced with the values from the analyzed spectrum,
								# such as to have a mass multiple times for multiple matched formulas
								
								mz <- c(msmsrecord$childFilt$mzFound,
										msmsrecord$childBad$mzFound,
										msmsrecord$childUnmatched$mzFound,
										msmsrecord$childRawLow[,mzcol], 
										msmsrecord$childRawSatellite[,mzcol])
								if(length(mz) > 0)
								{
									spectrum@mz <- mz
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
									spectrum@info <- msmsrecord$info
								}
								else
								{
									# mz and intensity are already there
									spectrum@satellite <- as.logical(rep(NA, spectrum@peaksCount))
									spectrum@low <- as.logical(rep(NA, spectrum@peaksCount))
									spectrum@rawOK <- as.logical(rep(FALSE, spectrum@peaksCount))
									spectrum@good <- as.logical(rep(FALSE, spectrum@peaksCount))
									spectrum@mzCalc <- as.numeric(rep(NA, spectrum@peaksCount))
									spectrum@formula <- as.character(rep(NA, spectrum@peaksCount))
									spectrum@dbe <- as.numeric(rep(NA, spectrum@peaksCount))
									spectrum@formulaCount <- as.integer(rep(NA, spectrum@peaksCount))
									spectrum@dppm <- as.numeric(rep(NA, spectrum@peaksCount))
									spectrum@dppmBest <- as.numeric(rep(NA, spectrum@peaksCount))
									
									
									
								}
								
								
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
							
							ok <- unlist(lapply(set@children, function(c) c@ok))
							if(all(ok))
								set@complete <- TRUE
							if(all(!ok))
								set@empty <- TRUE
							
							
							set@mode <- analyzedSpec$mode
							return(set)
				},
				spectra, analyzedSpecs)
	}


	return(as(spectra, "SimpleList"))
	
}

.updateObject.aggregated <- function(aggregatedSpecs, refilteredSpecs =  NULL)
{
  # Note: instead of rbind()ing peaksMatched and peaksUnmatched together from the start, they are treated separately
  # and rbind()ed at the end. This makes the process easier because everything can be substituted in by matching rownames
  # (note: for refilteredRcSpecs this doesn't work anymore, that's why it's treated separately)
	aggregatedSpecs$peaksUnmatched$dppmRc <- NA
	
	if(!is.null(aggregatedSpecs$peaksUnmatchedC))
	{
    # add noise columns
		aggregatedSpecs$peaksUnmatched <- addProperty(aggregatedSpecs$peaksUnmatched, "noise", "logical", TRUE)
		aggregatedSpecs$peaksMatched <- addProperty(aggregatedSpecs$peaksMatched, "noise", "logical", FALSE)
		# Set the peaks which are in peaksUnmatchedC to non-noise
		aggregatedSpecs$peaksUnmatched[match(rownames(aggregatedSpecs$peaksUnmatchedC), rownames(aggregatedSpecs$peaksUnmatched)),
				"noise"] <- FALSE	
		
		
		if(!is.null(aggregatedSpecs$peaksReanalyzed))
		{
      # to both matched and unmatched peaks, add the reanalysis result columns.
      # Into the unmatched peaks table, substitute the results (with the corresponding row names) from old-workspace peaksReanalyzed
      
			aggregatedSpecs$peaksUnmatched <- addProperty(aggregatedSpecs$peaksUnmatched, "reanalyzed.formula", "character")
			aggregatedSpecs$peaksUnmatched <- addProperty(aggregatedSpecs$peaksUnmatched, "reanalyzed.mzCalc", "numeric")
			aggregatedSpecs$peaksUnmatched <- addProperty(aggregatedSpecs$peaksUnmatched, "reanalyzed.dppm", "numeric")
			aggregatedSpecs$peaksUnmatched <- addProperty(aggregatedSpecs$peaksUnmatched, "reanalyzed.formulaCount", "numeric")
			aggregatedSpecs$peaksUnmatched <- addProperty(aggregatedSpecs$peaksUnmatched, "reanalyzed.dbe", "numeric")
			aggregatedSpecs$peaksUnmatched <- addProperty(aggregatedSpecs$peaksUnmatched, "matchedReanalysis", "logical", FALSE)
			
			aggregatedSpecs$peaksUnmatched[match(rownames(aggregatedSpecs$peaksReanalyzed), rownames(aggregatedSpecs$peaksUnmatched)),
					c("reanalyzed.formula","reanalyzed.mzCalc","reanalyzed.dppm","reanalyzed.formulaCount","reanalyzed.dbe")] <-
					aggregatedSpecs$peaksReanalyzed[,
					c("reanalyzed.formula","reanalyzed.mzCalc","reanalyzed.dppm","reanalyzed.formulaCount","reanalyzed.dbe")]
			aggregatedSpecs$peaksUnmatched$matchedReanalysis <- !is.na(aggregatedSpecs$peaksUnmatched$reanalyzed.dppm)
			
			aggregatedSpecs$peaksMatched <- addProperty(aggregatedSpecs$peaksMatched, "reanalyzed.formula", "character")
			aggregatedSpecs$peaksMatched <- addProperty(aggregatedSpecs$peaksMatched, "reanalyzed.mzCalc", "numeric")
			aggregatedSpecs$peaksMatched <- addProperty(aggregatedSpecs$peaksMatched, "reanalyzed.dppm", "numeric")
			aggregatedSpecs$peaksMatched <- addProperty(aggregatedSpecs$peaksMatched, "reanalyzed.formulaCount", "numeric")
			aggregatedSpecs$peaksMatched <- addProperty(aggregatedSpecs$peaksMatched, "reanalyzed.dbe", "numeric")
			aggregatedSpecs$peaksMatched <- addProperty(aggregatedSpecs$peaksMatched, "matchedReanalysis", "logical", NA)
		}
		
	}
	
	
	specs <- rbind(aggregatedSpecs$peaksMatched, aggregatedSpecs$peaksUnmatched)
	
	# index the spectra
	specs <- addProperty(specs, "index", "integer")
	if(nrow(specs) > 0)
		specs$index <- 1:nrow(specs)
	
  # remove the mz column if present - it formerly held the unrecalibrated mz value in recalibrated aggregated tables
	specs$mz <- NULL
	
  # rename the int column to intensity
	specs$intensity <- specs$int
	specs$int <- NULL
		
	return(specs)
}

.updateObject.refiltered <- function(w, specs, refilteredSpecs)
{
  # Filter peaks again, and redetermine the filterMultiplicity settings heuristically - 
  # it is much easier than substituting the info in from the refiltering table
  
  # Get what the lowest multiplicity was, and filter using this.
  multiplicityFilter <- min(refilteredSpecs$peaksOK$formulaMultiplicity, refilteredSpecs$peaksReanOK$formulaMultiplicity)
  
  
  w <- filterMultiplicity(
    w, archivename = NA, mode = NA, multiplicityFilter = multiplicityFilter)
  # aggregate again:
  w@aggregated <- aggregateSpectra(w@spectra, addIncomplete=TRUE)
  w <- processProblematicPeaks(w, "")
  
  return(w)
  
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