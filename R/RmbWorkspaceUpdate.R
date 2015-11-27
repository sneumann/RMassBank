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
					slot(w.parent, "recalibratedSpecs", check=FALSE) <- NULL
					slot(w.parent, "analyzedRcSpecs", check=FALSE) <- NULL
					slot(w.parent, "aggregatedRcSpecs", check=FALSE) <- NULL
					slot(w.parent, "reanalyzedRcSpecs", check=FALSE) <- NULL
					slot(w.parent, "refilteredRcSpecs", check=FALSE) <- NULL
					slot(w.old, "specs", check=FALSE) <- w.old@recalibratedSpecs
					slot(w.old, "analyzedSpecs", check=FALSE) <- w.old@analyzedRcSpecs
					w.parent.new <- updateObject(w.parent)
					w.new@parent <- w.parent.new
				}
				w.new@spectra <- .updateObject.spectra(w.old@specs, w.old@analyzedSpecs)
				if(7 %in% progress)
				{
					w.new@aggregated <- .updateObject.aggregated(w.old@reanalyzedRcSpecs)
				}
				else if(6 %in% progress)
				{
					w.new@aggregated <- .updateObject.aggregated(w.old@aggregatedRcSpecs)
				}
        else if(3 %in% progress)
        {
          w.new@aggregated <- .updateObject.aggregated(w.old@aggregatedSpecs)
        }
        
        if(8 %in% progress)
        {
          w.new@aggregated <- .updateObject.refiltered(w.new, w.new@aggregated, w.old@refilteredRcSpecs)
        }
        
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
				set@id <- as.character(as.integer(spec$id))
				set@formula <- spec$formula
				set@found <- as.logical(spec$foundOK)
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
  
  peaksFiltered <- filterPeaksMultiplicity(peaksMatched(specs),
                                           "formula", TRUE)
  
  
  peaksFilteredReanalysis <- 
    filterPeaksMultiplicity(specs[!is.na(specs$matchedReanalysis) & specs$matchedReanalysis,,drop=FALSE], "reanalyzed.formula", FALSE)
  
  
  
  specs <- addProperty(specs, "formulaMultiplicity", "numeric", 0)
  
  # Reorder the columns of the filtered peaks such that they match the columns
  # of the original aggregated table; such that the columns can be substituted in.
  
  peaksFiltered <- peaksFiltered[,colnames(specs)]
  peaksFilteredReanalysis <- peaksFilteredReanalysis[,colnames(specs)]
  
  # substitute into the parent dataframe
  specs[match(peaksFiltered$index,specs$index),] <- peaksFiltered
  specs[match(peaksFilteredReanalysis$index,specs$index),] <- peaksFilteredReanalysis
  
  multiplicityFilter <- min(refilteredSpecs$peaksOK$formulaMultiplicity, refilteredSpecs$peaksReanOK$formulaMultiplicity)
  
  specs <- addProperty(specs, "filterOK", "logical", FALSE)
  
  specs[specs$formulaMultiplicity > (multiplicityFilter - 1),"filterOK"] <- TRUE
  
  
  
  peaksReanOK <- specs[
    specs$filterOK & !is.na(specs$matchedReanalysis) & specs$matchedReanalysis,,drop=FALSE]
  
  # build M+H+ table
  mhsat <- lapply(w@spectra, function(s) c(s@id, findMz.formula(s@formula, "pH", 10, 0)$mzCenter))
  mhsat.df <- as.data.frame(do.call(rbind, mhsat))
  colnames(mhsat.df) <- c("cpdID", "mass")
  mhsat.df$mass <- as.numeric(as.character(mhsat.df$mass))
  
  # Kick the M+H+ satellites out of peaksReanOK:
  peaksReanOK$mzCenter <- mhsat.df[match(as.numeric(as.character(peaksReanOK$cpdID)), as.numeric(as.character(mhsat.df$cpdID))),"mass"]
    
  peaksReanBad <- peaksReanOK[
    !((peaksReanOK$mzFound < peaksReanOK$mzCenter - 1) |
        (peaksReanOK$mzFound > peaksReanOK$mzCenter + 1)),]
  specs[match(peaksReanBad$index, specs$index),"filterOK"] <- FALSE
  
  return(specs)
  
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