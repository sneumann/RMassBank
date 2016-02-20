#' Checks for isotopes in a \code{msmsWorkspace}
#' 
#' @param w A \code{msmsWorkspace} to work with.
#' @param mode \code{"pH", "pNa", "pM", "pNH4", "mH", "mM", "mFA"} for different ions 
#' 			([M+H]+, [M+Na]+, [M]+, [M+NH4]+, [M-H]-, [M]-, [M+FA]-).
#' @param intensity_cutoff The cutoff (as an absolute intensity value) under which isotopic peaks shouldn't be checked for or accepted as valid.
#'        Please note: The cutoff is not hard in the sense that it interacts with the intensity_precision parameter.
#' @param intensity_precision The difference that is accepted between the calculated and observed intensity of a possible isotopic peak. Further details down below.
#' @param conflict Either "isotopic"(Peak formulas are always chosen if they fit the requirements for an isotopic peak)
#' 		  or "strict"(Peaks are only marked as isotopic when there hasn't been a formula assigned before.)
#' @param isolationWindow Half of the width of the isolation window in Da 
#' @param evalMode Currently no function yet, but planned. Currently must be "complete"
#' @param plotSpectrum A boolean specifiying whether the spectrumshould be plotted
#' @param settings Options to be used for processing. Defaults to the options loaded via
#' 			\code{\link{loadRmbSettings}} et al. Refer to there for specific settings.
#' @details text describing parameter inputs in more detail.
#' \itemize{
#'  \item{\code{intensity_precision}}{This parameter determines how strict the intensity values should adhere to the calculated intensity in relation to the parent peak.
#'  Options for this parameter are \code{"none"}, where the intensity is irrelevant, \code{"low"}, which has an error margin of 70\% and \code{"high"}, where the error margin 
#'  is set to 35\%. The recommended setting is \code{"low"}, but can be changed to adjust to the intensity precision of the mass spectrometer.}
#  \item{\code{evalMode}}{This parameter sets what should be done after the isotopic check. The option "add" adds failpeaks if they are isotopic peaks of previously matched peaks.
#  "check" checks matched peaks with formulas for isotopes and removes them if no isotopic peaks have been found for all formulas. The formula is also
#  adjusted, if the one with matching isotopes doesn't have the lowest dppm. "complete" does both.}
#' }
#' @return The \code{msmsWorkspace} with annotated isolation peaks
#' @author Michael Stravs, Eawag <michael.stravs@@eawag.ch>
#' @author Erik Mueller, UFZ
#' @export
checkIsotopes <- function(w, mode = "pH", intensity_cutoff = 1000, intensity_precision = "low", conflict = "dppm", 
							isolationWindow = 2, evalMode = "complete", plotSpectrum = TRUE, settings = getOption("RMassBank")){

	# Load library and data
	requireNamespace("enviPat",quietly=TRUE)
	
	data("isotopes", package="enviPat")
	
	if(!(intensity_precision %in% c("none","low","high"))){
		stop('intensity_precision must be specified as either "none", "low" or "high"')
	}
	
	switch(intensity_precision, 
		none={
			precisionVal <- Inf
		},
		low={
			precisionVal <- 0.7
		},
		high={
			precisionVal <- 0.35
		}
	)
	
	# Load filtersettings
	filterSettings = settings$filterSettings
	
	# Assign formula additions according to code
	if(mode == "pH") {
		allowed_additions <- "H"
		mode.charge <- 1
	} else if(mode == "pNa") {
		allowed_additions <- "Na"
		mode.charge <- 1
	} else if(mode == "pM") {
		allowed_additions <- ""
		mode.charge <- 1
	} else if(mode == "mM") {
		allowed_additions <- ""
		mode.charge <- -1
	} else if(mode == "mH") {
		allowed_additions <- "H-1"
		mode.charge <- -1
	} else if(mode == "mFA") {
		allowed_additions <- "C2H3O2"
		mode.charge <- -1
	} else if(mode == "pNH4") {
		allowed_additions <- "NH4"
		mode.charge <- 1
	} else{
		stop("mode = \"", mode, "\" not defined")
	}
	
	if(!(evalMode %in% c("complete"))){
	    stop('evalMode must currently be specified as "complete"')
	}
	
	evalMode <- c("add")
	
	
	# "default" isotopes (i.e. those with the highest abundance)
	defIsotopes <- c("107Ag", "27Al", "40Ar", "75As", "197Au", "11B", "138Ba", "9Be", "209Bi",
	"79Br", "12C", "40Ca", "114Cd", "140Ce", "35Cl", "59Co", "52Cr", "133Cs", "63Cu", "164Dy",
	"166Er", "153Eu", "19F", "56Fe", "69Ga", "158Gd", "74Ge", "1H", "4He", "180Hf", "202Hg",
	"165Ho", "127I", "115In", "193Ir", "39K", "84Kr", "139La", "7Li", "175Lu", "24Mg", "55Mn",
	"98Mo", "14N", "23Na", "93Nb", "142Nd", "20Ne", "58Ni", "16O", "192Os", "31P", "231Pa",
	"208Pb", "106Pd", "141Pr", "195Pt", "85Rb", "187Re", "103Rh", "102Ru", "32S", "121Sb",
	"45Sc", "80Se", "28Si", "152Sm", "120Sn", "86Sr", "88Sr", "181Ta", "159Tb", "130Te",
	"232Th", "48Ti", "205Tl", "169Tm", "238U", "51V", "184W", "132Xe", "89Y", "174Yb",
	"64Zn", "90Zr")

	
	# Get the ppm limit from the settings
	ppmlimit <- 2.25 * filterSettings$ppmFine
	
	# Get the cutoff from the settings
	cut <- filterSettings$fineCut
	cut_ratio <- filterSettings$fineCutRatio

	# Extract matched and unmatched peaks from the aggregated peaks
	matchedPeaks <- peaksMatched(w)
	unmatchedPeaks <- peaksUnmatched(w)
	
	wEnv <- environment()
	
	# lapply over all runs
	lapply(w@spectra, function(spec){
		
		# Find parent formula and cpdID
		parent_formula <- add.formula(spec@formula, allowed_additions)
		id <- as.numeric(spec@id)
		specNum <- 0
		specEnv <- environment()

		# lapply over all extracted MS2 spectra
		lapply(spec@children, function(msmsdata){
			
			specEnv$specNum <- specEnv$specNum + 1
			
			# Extract currently relevant peaks
			currentMPeaks <- matchedPeaks[(matchedPeaks$cpdID == id) & (matchedPeaks$scan == msmsdata@acquisitionNum),,drop=FALSE]
			currentUPeaks <- unmatchedPeaks[(unmatchedPeaks$cpdID == id) & (unmatchedPeaks$scan == msmsdata@acquisitionNum),,drop=FALSE]
			
			if(nrow(currentMPeaks)){
				rownames(currentMPeaks) <- 1:nrow(currentMPeaks)
			} else {
				message(paste0("Compound ", id, " in spectrum #", specEnv$specNum," does not have matched peaks, so no isotopes can be found"))
				if(plotSpectrum){
					plot(currentMPeaks$mzFound, currentMPeaks$intensity,type="h", main=paste(id,findName(id)), col="black", xlab="m/z", ylab="intensity", lwd=3)
				}
				return(0)
			}
			
			if(nrow(currentUPeaks)){
				rownames(currentUPeaks) <- 1:nrow(currentUPeaks)
			}
			
			# Generate isotopic patterns of the matched peaks
			# sort out possible isotopic peaks according to
			# isolationwindow and intensity
			isoMPatterns <- apply(currentMPeaks, 1, .findPattern, defIsotopes = defIsotopes, intensity_cutoff = intensity_cutoff, 
									precisionVal = precisionVal, ppmlimit = ppmlimit, isolationWindow = isolationWindow)
			
			# Name the isopatterns
			names(isoMPatterns) <- currentMPeaks$formula
			
			
			# Which isotope patterns still have theoretical intensities above the cutoff?
			peaksToCheck <- which(as.logical(sapply(isoMPatterns,nrow)))
			
			if(!length(peaksToCheck)){
				message(paste0("The peaks of compound ", id, " in spectrum #", specEnv$specNum," are not intense enough to search for isotopic peaks"))
				if(plotSpectrum){
					plot(currentMPeaks$mzFound, currentMPeaks$intensity,type="h", main=paste(id,findName(id)), col="black", xlab="m/z", ylab="intensity", lwd=3)
				}
				return(0)
			}
			
			# Now, look for isotopic patterns in unmatched peaks with all specified parameters

			# Which peaks have no formula annotation within the dppm as of now?
			peaksNoAnnotation <- currentUPeaks[which(is.na(currentUPeaks$dppm) | abs(currentUPeaks$dppmBest) > ppmlimit/2.25),]
			
			# If there are any peaks without annotation:
			if(nrow(peaksNoAnnotation)){
				UPList <- .findMatchingPeaks(peaksNoAnnotation,isoMPatterns[peaksToCheck])
			} else{
				# If there are no peaks, fake an empty dataset
				UPList <- list(list(data.frame(dummy=character())),list(data.frame(dummy=character())))
			}
			
			# Generate matrix of peaks that should be added
			additionMatrix <- .peakReasoner(UPList, currentMPeaks)
			
            # If there is something in the matrix
            if(nrow(additionMatrix)){
                # Add all the peaks that could be annotated as isotopic
                wEnv$w@aggregated[additionMatrix$index,] <- additionMatrix
            }
            
            
			# If conflict is "strict", end here (And plot, maybe)
			if(conflict == "strict"){
				if(plotSpectrum){
					plot(currentMPeaks$mzFound, currentMPeaks$intensity,type="h", main=paste(id,findName(id)), col="black", xlab="m/z", ylab="intensity", lwd=3)
					if(nrow(additionMatrix)){
						points(additionMatrix$mzFound, additionMatrix$intensity,type="h", col="green", lwd=3)
					}
				}
				return(0)
			}
			
			
			# Now check the matched peaks for the patterns and put all peaks in a matrix
			MPList <- .findMatchingPeaks(currentMPeaks[currentMPeaks$filterOK,], isoMPatterns[peaksToCheck])
			
			# Generate matrix of peaks that should be corrected
			correctionMatrix <- .peakReasoner(MPList, currentMPeaks)
            
            # Where have isotopes been found/not found in matched peaks
			IsoPeaksMindex   <- which(sapply(MPList, function(x) any(sapply(x,nrow))))
			noIsoPeaksMindex <- which(!sapply(MPList, function(x) any(sapply(x,nrow))))
			
            # Where have isotopes been found/not found in unmatched peaks
            IsoPeaksUindex   <- which(sapply(UPList, function(x) any(sapply(x,nrow))))
			noIsoPeaksUindex <- which(!sapply(UPList, function(x) any(sapply(x,nrow))))
			
			if(conflict=="isotopic"){
			    if(nrow(correctionMatrix)){
			        
			        # Peaks that are changed but also seem to have isotopes themselves
			        confPeaksIndex <- which(correctionMatrix$index %in% currentMPeaks$index[peaksToCheck[IsoPeaksMindex]])
			        conflictedMatrix <- correctionMatrix[confPeaksIndex,,drop=FALSE]
			        
			        if(length(confPeaksIndex)){
			            correctionMatrix <- correctionMatrix[-confPeaksIndex,]
			        }
			        
			        
			        if(nrow(correctionMatrix)){
			            wEnv$w@aggregated[correctionMatrix$index,] <- correctionMatrix
			        }
			    }
			    
			    if(!("check" %in% evalMode)){
			        if(plotSpectrum){
			            plot(currentMPeaks$mzFound, currentMPeaks$intensity,type="h", main=paste(id,findName(id)), col="black", xlab="m/z", ylab="intensity", lwd=3)
			            if(nrow(additionMatrix)){
			                points(additionMatrix$mzFound, additionMatrix$intensity,type="h", col="green", lwd=3)
			            }
			            if(nrow(correctionMatrix)){
			                points(correctionMatrix$mzFound, correctionMatrix$intensity,type="h", col="yellow", lwd=3)
			            }
			            if(nrow(conflictedMatrix)){
			                points(conflictedMatrix$mzFound, conflictedMatrix$intensity,type="h", col="red", lwd=3)
			            }
			        }
			        return(0)
			    }
			}
			
			if("check" %in% evalMode){

			    if(plotSpectrum){
			        plot(currentMPeaks$mzFound, currentMPeaks$intensity,type="h", main=paste(id,findName(id)), xlim=c(min(currentMPeaks$mzFound), max(currentMPeaks$mzFound)), col="black", xlab="m/z", ylab="intensity", lwd=3)
			    }
			    
				indexList <- lapply(unique(currentMPeaks$mzFound[peaksToCheck]),function(mz){
					which(currentMPeaks$mzFound[peaksToCheck] == mz)
				})
				
				# Matched Peaks where no isotopes have been found
				hasNoIsoIndex <- unlist(sapply(indexList,function(indices){
					if(any(IsoPeaksMindex %in% indices) || any(IsoPeaksUindex %in% indices)){
						return(vector())
					} else{
						return(currentMPeaks$index[peaksToCheck[indices]])
					}
				}))
				
				# Matched Peaks which should have isotopes and are no isotopes themselves
				isNoIsoIndex <- currentMPeaks$index[peaksToCheck[
					which(!(currentMPeaks$mzFound[peaksToCheck] %in% correctionMatrix$mzFound))
				]]
				
				noIsoPeaksMatrix <- wEnv$w@aggregated[intersect(hasNoIsoIndex,isNoIsoIndex),]
				
				
				
			    if(nrow(noIsoPeaksMatrix)){
			        wEnv$w@aggregated[noIsoPeaksMatrix$index,]$good <- FALSE
			        wEnv$w@aggregated[noIsoPeaksMatrix$index,]$filterOK <- FALSE
			        wEnv$w@aggregated[noIsoPeaksMatrix$index,]$matchedReanalysis <- FALSE
			        if(plotSpectrum){
			            points(noIsoPeaksMatrix$mzFound, noIsoPeaksMatrix$intensity, type="h", col="orange", lwd=3)
			        }
			    }
			    
			    if("add" %in% evalMode && conflict=="isotopic" && plotSpectrum){
			        if(nrow(additionMatrix)){
			            points(additionMatrix$mzFound, additionMatrix$intensity,type="h", col="green", lwd=3)
			        }
			        if(nrow(correctionMatrix)){
			            points(correctionMatrix$mzFound, correctionMatrix$intensity,type="h", col="yellow", lwd=3)
						if(nrow(conflictedMatrix)){
							points(conflictedMatrix$mzFound, conflictedMatrix$intensity,type="h", col="red", lwd=3)
						}
			        }

			    }
			}
			return(0)
			
		})
	})
	return(w)
}

# Pattern finding and evaluation of these patterns inside the checkIsotopes function harmed readability and complicated debugging
# So modularize this function
.findPattern <- function(aggregateRow, defIsotopes, intensity_cutoff, precisionVal, ppmlimit, isolationWindow){
	
	# Find pattern and mass
	outp <- capture.output(pattern <- as.data.frame(enviPat::isopattern(aggregateRow["formula"], isotopes = isotopes)[[1]]))
    rm("outp")
	mass <- as.numeric(aggregateRow["mzCalc"])
	
	# Find index of monoisotopic molecule and
	# normalize abundance that monoisotopic molecule has always "1"
	mainIndex <- which.min(abs(pattern[,"m/z"]-mass))
	pattern[,"abundance"] <- round(pattern[,"abundance"] * (100/pattern[mainIndex,"abundance"]),3)
	
    # Add the formula of every isotope to every row in the pattern
	pattern <- .annotateFormulaToEnviPatTable(pattern, defIsotopes)
	
	# Sort pattern by abundance
	pattern <- pattern[order(pattern[,"abundance"],decreasing = T),][-mainIndex,]
	
	# Look up which patterns have their m/z inside the isolation window +- 0.5
	# and delete all others 
	keepMzIndex <- which(pattern[,"m/z"] < mass + isolationWindow + 0.2 & pattern[,"m/z"] > mass - isolationWindow - 0.2)
	pattern <- pattern[keepMzIndex,,drop=FALSE]
	if(nrow(pattern) == 0){
		return(pattern)
	}
	
	# Calculate the expected intensities according to the abundance
	# See which expected isotope peaks have an intensity higher than the cutoff
	
	# Find the absolute intensity ranges
	intensities <- apply(pattern, 1, function(patternRow){
		as.numeric(aggregateRow["intensity"]) * as.numeric(patternRow["abundance"])/100
	})
	
	roundedInt <- round(intensities,digits=-2)
	keepIntIndex <- which((roundedInt + roundedInt * precisionVal) >= intensity_cutoff)
	pattern <- pattern[keepIntIndex,,drop=FALSE]
	if(nrow(pattern)){
		pattern$minintensity <- roundedInt[keepIntIndex] - roundedInt[keepIntIndex] * precisionVal
		pattern$maxintensity <- roundedInt[keepIntIndex] + roundedInt[keepIntIndex] * precisionVal
        pattern$expectedintensity <- roundedInt[keepIntIndex]
		
		# Calculate the expected mz range of the isotope peak
		mzCols <- t(sapply(pattern[,"m/z"], ppm, dppm=ppmlimit/2.25,l=T))
		colnames(mzCols) <- c("maxMZ","minMZ")
		pattern <- cbind(pattern,mzCols)
	}
	
	return(pattern)
}

.findMatchingPeaks <- function(aggregatedPeakMatrix, isoPatterns){
	# Only unique mzs (peaks can be annotated multiple times and will therefore be included multiple times)
	# It doesn't matter which one, so take the first
	aggregatedPeakMatrix <- aggregatedPeakMatrix[sapply(unique(aggregatedPeakMatrix$mzFound), function(mz){
		which(aggregatedPeakMatrix$mzFound == mz)[1]
	}),]
	
	# Iterate over all supplied patterns
	lapply(isoPatterns, function(pattern){
		x <- list()
		
		apply(pattern, 1, function(patternRow){
			# Find peaks that fit the specified intensities and m/z
			pIndex <- which(aggregatedPeakMatrix[,"mzFound"] < as.numeric(patternRow["maxMZ"])
									& aggregatedPeakMatrix[,"mzFound"] > as.numeric(patternRow["minMZ"])
									& aggregatedPeakMatrix[,"intensity"] < as.numeric(patternRow["maxintensity"])
									& aggregatedPeakMatrix[,"intensity"] > as.numeric(patternRow["minintensity"])
			)
			
			# Note these Peaks
			candidates <- aggregatedPeakMatrix[pIndex,]
			# If there are any: Change parameters in the aggregated matrix
			if(nrow(candidates)){
				candidates$dppm <- (candidates$mzFound /  as.numeric(patternRow["m/z"]) - 1) * 1e6
				candidates$mzCalc <-  as.numeric(patternRow["m/z"])
				candidates$formula <-  patternRow["formula"]
				candidates$matchedReanalysis <- NA
				candidates$filterOK <- TRUE
				candidates$good <- TRUE
				candidates$dppmBest <- candidates$dppm
				candidates$formulaCount <- 1
			}
			return(candidates)
		})
	})
}

.peakReasoner <- function(adjustmentList, aggregatedPeakMatrix){
				
	adjustmentMatrix <- do.call(rbind, unlist(adjustmentList,recursive=FALSE))
	
	if(!as.logical(nrow(adjustmentMatrix))){
		return(adjustmentMatrix)
	}
	
	# Peaks can turn up multiple times, take the one that actually fits the "base" formula,
	# i.e. that of the peak that is actually written into the record
	
	# For each possible parent peak: Find out which peaks have possible children
	# and how many, group these peaks (by adding a column to the adjustmentMatrix)
	# Note: Grouping happens here, because we also need to split the peaks by mz
	summaryList <- lapply(adjustmentList, function(x) sapply(x,nrow))
	groupLengths <- unlist(lapply(summaryList, function(entry){
		if(all(entry == 0)) return(vector()) else return(length(which(entry > 0)))
	}))
	adjustmentMatrix$group <- unlist(lapply(1:length(groupLengths),function(i) rep(i,groupLengths[i])))
	
	# Technically, you only need to look at one peak in each group and can adjust the rest alongside
	# Since we only want to check if the parent formula is currently passed through
	groupAdjustmentMatrix <- adjustmentMatrix[sapply(1:length(groupLengths),function(i) which(adjustmentMatrix$group == i)[1]),]
	
	# Extract the parent formulas of the isotopic formula groups
	parentFormulas <- names(groupLengths)
	
	# Find out which of the possible found isotopes formulas are competing (by checking which have the same mz)
	indexList <- lapply(unique(groupAdjustmentMatrix$mzFound),function(mz){
		which(groupAdjustmentMatrix$mzFound == mz)
	})
	
	# For every set of indices, find the parent formula in the aggregated matrix
	return(do.call(rbind, lapply(indexList, function(indices){
		# All peaks which fit the parent formulas
		cutAggregateMatrix <- aggregatedPeakMatrix[which(aggregatedPeakMatrix$formula %in% parentFormulas[indices]),,drop=FALSE]
		correctParentIndex <- which(cutAggregateMatrix$filterOK == TRUE)
		if(length(correctParentIndex)){
			# If there is a Peak that will be passed and fits one of the parent formulas, take that Peak and all belonging to that group
			# Also, drop the group column
			correctGroup <- groupAdjustmentMatrix$group[indices[correctParentIndex]]
			returnMatrix <- adjustmentMatrix[which(adjustmentMatrix$group == correctGroup),-25,drop=FALSE]
			return(returnMatrix)
		} else{
			# Else calculate the Peak-isotope relationship with the lowest average dppm
			# And correct the Peak that currently has filterOK to not have it
			correctParentIndex <- which.min((abs(cutAggregateMatrix$dppm) + abs(adjustmentMatrix[indices,]$dppm))/2)
			correctionLine <- aggregatedPeakMatrix[which(aggregatedPeakMatrix$mzFound == groupAdjustmentMatrix$mzFound[indices[1]] & aggregatedPeakMatrix$filterOK),,drop=FALSE]
			correctionLine$filterOK <- FALSE
			correctGroup <- groupAdjustmentMatrix$group[indices[correctParentIndex]]
			returnMatrix <- adjustmentMatrix[which(adjustmentMatrix$group == correctGroup),-25,drop=FALSE]
			return(rbind(
				returnMatrix,
				correctionLine
			))
		}
	})))
}


# Further modularization of formula addition
# Is essentially one task, so this makes it more understandable
.annotateFormulaToEnviPatTable <- function(pattern, defIsotopes){
    # Find all isotope atom names (only in colnames of patterns, sadly)
	isoCols <- colnames(pattern)[3:ncol(pattern)]
	
	# Order the formulas to be in Hill system:
	# First the Cs, then the Hs, then lexicographical
	# First: Split isotope names so that there only are atoms
	# Find position of C and H
	splitNames <- gsub('^([0-9]+)(.+)$', '\\2', isoCols)
	CHPos <- c(which(splitNames == "C"),which(splitNames == "H"))
	
	# Account for special case: No Cs or Hs
	if(length(CHPos)){
		# If there are Cs and Hs, overwrite the positions so they always get sorted to the front
		splitNames[CHPos] <- ""
	}
	
	# Default isotopes are always first in the colnames, so no need to order internally between isotopes
	atomOrder <- order(splitNames)
	
	# If there are Cs and Hs, the order must be preserved
	if(length(CHPos)){
		atomOrder[1:length(CHPos)] <- CHPos
	}
	# else it is already ordered
	
	# Mark new names for formula creation.
	newnames <- unname(sapply(isoCols, function(currentCol){
		if(currentCol %in% defIsotopes){
			# "Default" isotope (no isotope mass in formula, e.g. "C")
			return(gsub('^(.{0})([0-9]+)(.+)$', '\\3', currentCol))
		}else{
			# Other isotopes (isotope mass in formula, e.g. "[13]C")
			return(gsub('^(.{0})([0-9]+)(.+)$', '\\1[\\2]\\3', currentCol))
		}
	}))
	
	# Generate the formula for every pattern
	pattern$formula <- apply(pattern,1,function(p){
		paste0(sapply(atomOrder+2,function(isoIndex){
			if(p[isoIndex] == 0){ 
				return("")
			}
			if(p[isoIndex] == 1){
				return(newnames[isoIndex-2])
			}
			return(paste0(newnames[isoIndex-2],p[isoIndex]))
		}),collapse="")
	})
    
    return(pattern)
}