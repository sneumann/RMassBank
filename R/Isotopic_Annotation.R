#' Checks for isotopes in a \code{msmsWorkspace}
#' 
#' @param w A \code{msmsWorkspace} to work with.
#' @param mode \code{"pH", "pNa", "pM", "pNH4", "mH", "mM", "mFA"} for different ions 
#' 			([M+H]+, [M+Na]+, [M]+, [M+NH4]+, [M-H]-, [M]-, [M+FA]-).
#' @param intensity_cutoff The cutoff (as an absolute intensity value) under which isotopic peaks shouldn't be checked for or accepted as valid.
#' @param intensity_precision The difference that is accepted between the calculated and observed intensity of a possible isotopic peak. Further details down below.
#' @param conflict Either "isotopic"(Peak formulas are always chosen if they fit the requirements for an isotopic peak)
#' 		  or "strict"(Peaks are only marked as isotopic when there hasn't been a formula assigned before.)
#' @param isolationWindow The width of the isolation window in Da 
#' @param evalMode Currently no function yet, but planned
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
checkIsotopes <- function(w, mode = "pH", intensity_cutoff = 5000, intensity_precision = "low", conflict = "dppm", 
							isolationWindow = 4, evalMode = "complete", plotSpectrum = TRUE, settings = getOption("RMassBank")){

	# Load library and data
	require(enviPat)
	data("isotopes")
	
	
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
	
	evalMode <- c("add","check")
	
	
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
		

		# lapply over all extracted MS2 spectra
		lapply(spec@children, function(msmsdata){
			
			# Extract currently relevant peaks
			currentMPeaks <- matchedPeaks[(matchedPeaks$cpdID == id) & (matchedPeaks$scan == msmsdata@acquisitionNum),,drop=FALSE]
			currentUPeaks <- unmatchedPeaks[(unmatchedPeaks$cpdID == id) & (unmatchedPeaks$scan == msmsdata@acquisitionNum),,drop=FALSE]
			
			rownames(currentMPeaks) <- 1:nrow(currentMPeaks)
			rownames(currentUPeaks) <- 1:nrow(currentUPeaks)
			
			# Generate isotopic patterns of the matched peaks
			# Sort pattern by abundance
			isoMPatterns <- lapply(currentMPeaks$formula, function(formula){
				# Find pattern
				pattern <- as.data.frame(isopattern(formula, isotopes = isotopes)[[1]])
				mass <- findMz.formula(formula,"")$mzCenter
				
				# Find index of nonisotopic molecule and
				# normalize abundance that nonisotopic molecule has always "1"
				mainIndex <- which.min(abs(pattern[,"m/z"]-mass))
				pattern[,"abundance"] <- round(pattern[,"abundance"] * (100/pattern[mainIndex,"abundance"]),3)
				
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
					# If there are Cs and Hs, overwrite the positions so the always get sorted to the front
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
				pattern <- pattern[order(pattern[,"abundance"],decreasing = T),][-mainIndex,]
			})
			
			names(isoMPatterns) <- currentMPeaks$formula
			
			
			# Order patterns by abundance and make sure that the normal abundance is at 100%
			
			# Calculate the expected intensities according to the abundance
			# See which expected isotope peaks have an intensity higher than the cutoff
			
			for(foundAnnotation in 1:nrow(currentMPeaks)){
				intensities <- vector()
				
				for(pattern in 1:nrow(isoMPatterns[[foundAnnotation]])){
					intensities[pattern]  <- currentMPeaks$int[foundAnnotation] * isoMPatterns[[foundAnnotation]]$abundance[pattern]/100
				}
				
				roundedInt <- round(intensities,digits=-2)
				keepIndex <- which(roundedInt >= intensity_cutoff)
				isoMPatterns[[foundAnnotation]] <- isoMPatterns[[foundAnnotation]][keepIndex,,drop=FALSE]
				if(nrow(isoMPatterns[[foundAnnotation]])){
					isoMPatterns[[foundAnnotation]]$minintensity <- roundedInt[keepIndex] - roundedInt[keepIndex] * precisionVal
					isoMPatterns[[foundAnnotation]]$maxintensity <- roundedInt[keepIndex] + roundedInt[keepIndex] * precisionVal
					# Calculate the expected mz range of the isotope peak
					isoMPatterns[[foundAnnotation]] <- cbind(isoMPatterns[[foundAnnotation]],t(sapply(isoMPatterns[[foundAnnotation]][,"m/z"], ppm, dppm=ppmlimit,l=T)))
				}
			}
			
			# Which isotope patterns still have theoretical intensities above the cutoff?
			peaksToCheck <- which(as.logical(sapply(isoMPatterns,nrow)))
			
			# Now, look for isotopic patterns in unmatched peaks with all specified parameters
			if("add" %in% evalMode){
				# Which peaks have no correct formula annotation as of now?
				currentUPeaks <- currentUPeaks[which(is.na(currentUPeaks$dppm) | abs(currentUPeaks$dppmBest) > ppmlimit/2.25),]
				
				# If there are any peaks without annotation:
				if(nrow(currentUPeaks)){
					j <- 1
					UPList <- list()
					# Iterate over all patterns that have relevant intensities
					for(i in peaksToCheck){
						UPList[[j]] <- list()
						# Iterate over every isotopic permutation that is still in the pattern
						for(k in 1:nrow(isoMPatterns[[i]])){
							# Find peaks that fit the specified intensities and m/z
							pIndex <- which(currentUPeaks[,"mzFound"] < isoMPatterns[[i]][k,"1"] & currentUPeaks[,"mzFound"] > isoMPatterns[[i]][k,"2"]
													& currentUPeaks[,"intensity"] < isoMPatterns[[i]][k,"maxintensity"] & currentUPeaks[,"intensity"] > isoMPatterns[[i]][k,"minintensity"]
												)
							# Note these Peaks
							UPList[[j]][[k]] <- currentUPeaks[pIndex,]
							# If there are any: Change parameters in the aggregated matrix
							if(nrow(UPList[[j]][[k]])){ 
								UPList[[j]][[k]]$dppm <- (currentUPeaks[pIndex,]$mzFound / isoMPatterns[[i]][k,"m/z"] - 1) * 1e6
								UPList[[j]][[k]]$mzCalc <- isoMPatterns[[i]][k,"m/z"]
								UPList[[j]][[k]]$formula <- isoMPatterns[[i]][k,"formula"]
								UPList[[j]][[k]]$matchedReanalysis <- NA
								UPList[[j]][[k]]$filterOK <- TRUE
								UPList[[j]][[k]]$good <- TRUE
								UPList[[j]][[k]]$dppmBest <- UPList[[j]][[k]]$dppm
								UPList[[j]][[k]]$formulaCount <- 1
							}
						}
						j <- j + 1
					}
				} else{
					# If there are no peaks, fake an empty dataset
					UPList <- list(list(data.frame()),list(data.frame()))
				}
				
				
				additionMatrix <- do.call(rbind, unlist(UPList,recursive=FALSE))
				if(nrow(additionMatrix)){
					# Add all the peaks that could be annotated as isotopic
					wEnv$w@aggregated[additionMatrix$index,] <- additionMatrix
				}
				
				# If conflict is "strict", end here
				if(conflict == "strict"){
					if(plotSpectrum){
						plot(currentMPeaks$mzFound, currentMPeaks$intensity,type="h", main=paste(id,findName(id)), col="black", xlab="m/z", ylab="intensity", lwd=3)
						if(nrow(additionMatrix)){
							points(additionMatrix$mzFound, additionMatrix$intensity,type="h", col="green", lwd=3)
						}
					}
					return(0)
				}
			}
			
			
			# Now check the matched peaks for the patterns
			j <- 1
			MPList <- list()
			for(i in peaksToCheck){
				MPList[[j]] <- list()
				for(k in 1:nrow(isoMPatterns[[i]])){
					pIndex <- which(currentMPeaks[,"mzFound"] < isoMPatterns[[i]][k,"1"] & currentMPeaks[,"mzFound"] > isoMPatterns[[i]][k,"2"]
											& currentMPeaks[,"intensity"] < isoMPatterns[[i]][k,"maxintensity"] & currentMPeaks[,"intensity"] > isoMPatterns[[i]][k,"minintensity"]
											& currentMPeaks[,"filterOK"])
					MPList[[j]][[k]] <- currentMPeaks[pIndex,]

					if(nrow(MPList[[j]][[k]])){ 
						MPList[[j]][[k]]$dppm <- (currentMPeaks[pIndex,]$mzFound / isoMPatterns[[i]][k,"m/z"] - 1) * 1e6
						MPList[[j]][[k]]$mzCalc <- isoMPatterns[[i]][k,"m/z"]
						MPList[[j]][[k]]$formula <- isoMPatterns[[i]][k,"formula"]
						MPList[[j]][[k]]$matchedReanalysis <- NA
						MPList[[j]][[k]]$filterOK <- TRUE
						MPList[[j]][[k]]$good <- TRUE
						MPList[[j]][[k]]$dppmBest <- MPList[[j]][[k]]$dppm
						MPList[[j]][[k]]$formulaCount <- 1
					}
				}
				j <- j + 1
			}
			
			# Peakindices of peaksToCheck where no isotopes have been found in Unmatched Peaks
			noIsoPeaksUindex <- which(!sapply(UPList, function(x) any(sapply(x,nrow))))
			
			# Peakindices of peaksToCheck where isotopes have been found in Unmatched Peaks
			IsoPeaksUindex <- setdiff(1:length(peaksToCheck),noIsoPeaksUindex)
			
			# Peakindices of peaksToCheck where no isotopes have been found in Matched Peaks
			noIsoPeaksMindex <- which(!sapply(MPList, function(x) any(sapply(x,nrow))))
			
			# Peakindices of peaksToCheck where isotopes have been found in Matched Peaks
			IsoPeaksMindex <- setdiff(1:length(peaksToCheck),noIsoPeaksMindex)
			
			if("add" %in% evalMode && conflict=="isotopic"){
			    correctionMatrix <- do.call(rbind, unlist(MPList,recursive=FALSE))
			    correctionMatrix <- correctionMatrix[order(abs(correctionMatrix$dppm)),]
			    for(ind in unique(correctionMatrix$index)){
			        if(length(which(correctionMatrix$index==ind))>1){
			            correctionMatrix <- correctionMatrix[-which(correctionMatrix$index==ind)[-1],]
			        }
			    }
			    
			    if(nrow(correctionMatrix)){
			        
			        # Peaks that are changed but also seem to have isotopes themselves
			        confPeaksIndex <- which(correctionMatrix$index %in% currentMPeaks[peaksToCheck[IsoPeaksMindex],"index"])
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
			    # Matched Peaks where no isotopes have been found
			    noIsoPeaksMatrix <- currentMPeaks[peaksToCheck[noIsoPeaksMindex[which(noIsoPeaksMindex %in% noIsoPeaksUindex)]],]
			    
			    # Matched Peaks where no isotopes have been found and which haven't been marked as isotopic themselves
			    noIsoPeaksMatrix <- noIsoPeaksMatrix[which(!grepl("[",wEnv$w@aggregated$formula[noIsoPeaksMatrix$index],fixed=TRUE)),]
			    if(plotSpectrum){
			        plot(currentMPeaks$mzFound, currentMPeaks$intensity,type="h", main=paste(id,findName(id)), col="black", xlab="m/z", ylab="intensity", lwd=3)
			    }
			    
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
			        }
			        if(nrow(conflictedMatrix)){
			            points(conflictedMatrix$mzFound, conflictedMatrix$intensity,type="h", col="red", lwd=3)
			        }
			    }
			}
			return(0)
			
		})
	})
	return(w)
}
