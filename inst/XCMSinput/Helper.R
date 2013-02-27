##Helperscript to test a few things.
##This should work on any PC

library(RMassBank)
############################
############################
##TRYING TO GET THIS TO WORK
############################
############################

  analyzeTandemShot <- function(shot_mat)
  {
    shot <- as.data.frame(shot_mat)
    shot_orig <- shot
    # Filter out low intensity peaks:
    shot_lo <- shot[(shot$int < cut) | (shot$int < max(shot$int)*cut_ratio),]
    shot <- shot[(shot$int >= cut) & (shot$int > max(shot$int) * cut_ratio),]
    shot_full <- shot
    
    # Is there still anything left?
    if(nrow(shot)==0)
      return(list(specOK=FALSE))
    
    # Filter out satellite peaks:
    shot <- filterPeakSatellites(shot)
    shot_satellite_n <- setdiff(row.names(shot_full), row.names(shot))
    shot_satellite <- shot_full[shot_satellite_n,]

    # Is there still anything left?
    if(nrow(shot)==0)
      return(list(specOK=FALSE))
    
    if(max(shot$int) < 1e4)
      return(list(specOK=FALSE))
    # Crop to 4 digits (necessary because of the recalibrated values)
    shot[,mzColname] <- round(shot[,mzColname], 5)
    
	
	# here follows the Rcdk analysis
	#------------------------------------
	parentPeaks <- data.frame(mzFound=msmsPeaks$mz$mzCenter, 
			formula=msmsPeaks$formula,
			dppm=0,
			x1=0,x2=0,x3=0)
	
	# define the adduct additions
	if(mode == "pH")
	{
		allowed_additions <- "H"
		mode.charge <- 1
	}
	if(mode == "pNa")
	{
		allowed_additions <- "Na"
		mode.charge <- 1
	}
	if(mode == "pM")
	{
		allowed_additions <- ""
		mode.charge <- 1
	}
	if(mode == "mH")
	{
		allowed_additions <- "H-1"
		mode.charge <- -1
	}
	if(mode == "mFA")
	{
		allowed_additions <- "C2H3O2"
		mode.charge <- -1
	}
	
	# the ppm range is two-sided here.
	# The range is slightly expanded because dppm calculation of
	# generate.formula starts from empirical mass, but dppm cal-
	# culation of the evaluation starts from theoretical mass.
	# So we don't miss the points on 'the border'.
	
	if(run=="preliminary")
		ppmlimit <- 40
	else
		ppmlimit <- 15
	parent_formula <- add.formula(msmsPeaks$formula, allowed_additions)
	dbe_parent <- dbe(parent_formula)
	# check whether the formula is valid, i.e. has no negative or zero element numbers.
	#print(parent_formula)
	if(!is.valid.formula(parent_formula))
		return(list(specOK=FALSE))
	limits <- to.limits.rcdk(parent_formula)
	
	
	peakmatrix <- lapply(shot[,mzColname], function(mass) {
				peakformula <- tryCatch(generate.formula(mass, ppm(mass, ppmlimit, p=TRUE), 
								limits, charge=mode.charge), error=function(e) NA)
				#peakformula <- tryCatch( 
				#  generate.formula(mass, 
				#                   ppm(mass, ppmlimit, p=TRUE),
				#                   limits, charge=1),
				#error= function(e) list())
				if(!is.list(peakformula))
					return(t(c(mzFound=as.numeric(as.character(mass)),
											formula=NA, mzCalc=NA)))
				else
				{
					return(t(sapply(peakformula, function(f)
											{
												c(mzFound=mass,
														formula=f@string, 
														mzCalc=f@mass)
											})))
				}
			})
	
	childPeaks <- as.data.frame(do.call(rbind, peakmatrix))
	childPeaks$mzFound <- as.numeric(as.character(childPeaks$mzFound))
	childPeaks$formula <- as.character(childPeaks$formula)
	childPeaks$mzCalc <- as.numeric(as.character(childPeaks$mzCalc))
	childPeaks$dppm <- (childPeaks$mzFound / childPeaks$mzCalc - 1) * 1e6
	# delete the NA data out again, because MolgenMsMs doesn't have them
	# here and they will be re-added later
	# (this is just left like this for "historical" reasons)
	childPeaks <- childPeaks[!is.na(childPeaks$formula),]
	# check if a peak was recognized (here for the first time,
	# otherwise the next command would fail)
	if(nrow(childPeaks)==0)
		return(list(specOK=FALSE))
	
	# now apply the rule-based filters to get rid of total junk:
	# dbe >= -0.5, dbe excess over mother cpd < 3
	childPeaks$dbe <- unlist(lapply(childPeaks$formula, dbe))
	#iff_rcdk_pM_eln$maxvalence <- unlist(lapply(diff_rcdk_pM_eln$formula.rcdk, maxvalence))
	childPeaks <- childPeaks[childPeaks$dbe >= -0.5,] # & dbe < dbe_parent + 3)
	
	# check if a peak was recognized
	if(nrow(childPeaks)==0)
		return(list(specOK=FALSE))
	
	# trim mz to 5 digits
	shot[,mzColname] <- round(shot[,mzColname], 5)
    
    childPeaksInt <- merge(childPeaks, shot, by.x = "mzFound", by.y = mzColname, all.x = TRUE, all.y = FALSE )
    # find the best ppm value
    bestPpm <- aggregate(childPeaksInt$dppm, list(childPeaksInt$mzFound),
                         function(dppm) dppm[[which.min(abs(dppm))]])
    colnames(bestPpm) <- c("mzFound", "dppmBest")
    childPeaksInt <- merge(childPeaksInt, bestPpm, by="mzFound", all.x=TRUE)
    # count formulas found per mass
    countFormulasTab <- xtabs( ~formula + mzFound, data=childPeaksInt)
    countFormulas <- colSums(countFormulasTab)
    childPeaksInt$formulaCount <- countFormulas[as.character(childPeaksInt$mzFound)]
    # filter results
    childPeaksFilt <- filterLowaccResults(childPeaksInt, filterMode)
    childPeaksGood <- childPeaksFilt[["TRUE"]]
    childPeaksBad <- childPeaksFilt[["FALSE"]]
    # count formulas within new limits
    # (the results of the "old" count stay in childPeaksInt and are returned
    # in $childPeaks)
    if(!is.null(childPeaksGood))
    {
      countFormulasTab <- xtabs( ~formula + mzFound, data=childPeaksGood)
      countFormulas <- colSums(countFormulasTab)
      childPeaksGood$formulaCount <- countFormulas[as.character(childPeaksGood$mzFound)]
    }
    
    childPeaksUnmatched <- merge(childPeaks, shot, by.x = "mzFound", by.y = mzColname, 
                                 all.x = TRUE, all.y = TRUE )
    childPeaksUnmatched$dppmBest <- NA
    childPeaksUnmatched$formulaCount <- 0
    childPeaksUnmatched$good <- FALSE
    childPeaksUnmatched <- childPeaksUnmatched[is.na(childPeaksUnmatched$mzCalc),]
    
    # return list:
    rl <- list(
      specOK = !is.null(childPeaksGood),
      parent = parentPeaks,
      childFilt = childPeaksGood,
      childRaw=shot_orig
      )
    # if "detail" is set to TRUE, return more detailed results including
    # all the deleted peaks and the stages when they were culled
    if(detail)
    {
      rl$childRawLow = shot_lo
      rl$childRawSatellite = shot_satellite
      rl$childRawOK = shot
      rl$child =childPeaksInt
      rl$childBad = childPeaksBad
      rl$childUnmatched = childPeaksUnmatched
    }
    return(rl)
  }
  
  
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 
  
  
  
  
  
  
 XCMStoRMB <- function(msmsXCMSspecs, cpdID){
	ret <- list()
	ret$foundOK <- 1
	
	##Write nothing in the parents
	ret$parentscan <- 1
	ret$parentHeader <- matrix(0, ncol = 20, nrow = 1)
	rownames(ret$parentHeader) <-1
	colnames(ret$parentHeader) <- c("seqNum", "acquisitionNum", "msLevel", "peaksCount", "totIonCurrent", "retentionTime", "basepeakMZ", 
									"basePeakIntensity", "collisionEnergy", "ionisationEnergy", "lowMZ", "highMZ", "precursorScanNum",
									"precursorMZ", "precursorCharge", "precursorIntensity", "mergedScan", "mergedResultScanNum", 
									"mergedResultStartScanNum", "mergedResultEndScanNum")
	ret$parentHeader[1,1:3] <- 1
	ret$parentHeader[1,4:20] <- 0
	ret$parentHeader <- as.data.frame(ret$parentHeader)
	
	##Write the peaks into the childscans
	ret$childScans <- 2
	ret$childHeader <- matrix(0, ncol = 20, nrow = 1)
	rownames(ret$childHeader) <- 2
	colnames(ret$childHeader) <- c("seqNum", "acquisitionNum", "msLevel", "peaksCount", "totIonCurrent", "retentionTime", "basepeakMZ", 
									"basePeakIntensity", "collisionEnergy", "ionisationEnergy", "lowMZ", "highMZ", "precursorScanNum",
									"precursorMZ", "precursorCharge", "precursorIntensity", "mergedScan", "mergedResultScanNum", 
									"mergedResultStartScanNum", "mergedResultEndScanNum")
	ret$childHeader[1,1:2] <- ret$childScans
	ret$childHeader[1,3] <- 2
	ret$childHeader[1,4] <- length(msmsXCMSspecs[,1])
	ret$childHeader[1,5] <- 0 ##Does this matter?
	ret$childHeader[1,6] <- median(msmsXCMSspecs[,4])
	ret$childHeader[1,7] <- msmsXCMSspecs[which.max(msmsXCMSspecs[,7]),1]
	ret$childHeader[1,8] <- max(msmsXCMSspecs[,7])
	ret$childHeader[1,9] <- 0 ##Does this matter?
	ret$childHeader[1,10] <- 0 ##Does this matter?
	ret$childHeader[1,11] <- min(msmsXCMSspecs[,1])
	ret$childHeader[1,12] <- max(msmsXCMSspecs[,1])
	ret$childHeader[1,13] <- 1 
	ret$childHeader[1,14] <- findMz(cpdID)[[3]]
	ret$childHeader[1,15] <- 1 ##Will be changed for different charges
	ret$childHeader[1,16] <- 0 ##There sadly isnt any precursor intensity to find in the msms-scans. WorkarmsmsXCMS@files[1]ound?
	ret$childHeader[1,17:20] <- 0 ##Will be changed if merge is wanted
	ret$childHeader <- as.data.frame(ret$childHeader)
	ret$parentPeak <- matrix(nrow = 1, ncol = 2)
	colnames(ret$parentPeak) <- c("mz","int")
	ret$peaks <- list()
	ret$peaks[[1]] <- matrix(nrow = length(msmsXCMSspecs[,1]), ncol = 2)
	colnames(ret$peaks[[1]]) <- c("mz","int")
	ret$peaks[[1]][,1] <- msmsXCMSspecs[,1]
	ret$peaks[[1]][,2] <- msmsXCMSspecs[,7]
	ret$mz <- findMz(cpdID)
	ret$formula <- findFormula(cpdID)
	return(ret)
}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
##This uses an already created template for the settings
##Change path or settings at own risk
loadRmbSettings(system.file("XCMSinput/mysettings.ini",package="RMassBank"))

##Use the xcms-CAMERA-Peakpicker with the ChelidonineMSn-mzData
##1) Extract filepath
##2) Load Compoundlist
##3) Get the spec into specsXCMS
msmsXCMS <- newMsmsWorkspace()
filesXCMS <- system.file("XCMSinput/Chelidonine_666_pos.mzData",package="RMassBank") 
msmsXCMS@files <- filesXCMS
loadList(system.file("XCMSinput/Chelidonine.csv",package="RMassBank"))

##
## More example data
##
## http://www.casmi-contest.org/challenges-cat1-2.shtml
## Challenge3 MSMSneg10_Challenge3 MSMSneg20_Challenge3 MSMSneg30_Challenge3 MSMSneg40_Challenge3
##    4 individuelle Rohdatenfiles
##    Modus "mH" !
## http://www.casmi-contest.org/solutions-cat1-2.shtml
## Glucolesquerellin

##
## Overwrite parameters in spectraList
## 
getOption("RMassBank")$spectraList

## spectraList:
##  # First scan: CID 20
## - mode: CID
##   ces: 20
##   ce: 20
##   res: ca. 7500

xr <- system.file("microtofq/MSMSpos20_6.mzML", package="msdata", includeMsn=TRUE)
str(xr@msnCollisionEnergy)

msmsXCMSspecs <- findMsMsHRperxcms()
msmsXCMS@specs[[1]] <- XCMStoRMB(msmsXCMSspecs,666)
names(msmsXCMS@specs) <- findName(666)


mode = "pH"
#shots <- lapply(msmsXCMS@specs[[1]]$peaks, analyzeTandemShot)
msmsXCMS@analyzedSpecs <- lapply(msmsXCMS@specs, function(spec) {
				  s <- analyzeMsMs(spec, mode=mode, detail=TRUE, run="preliminary" )
				  return(s)
			  })


#########################################
#########################################
#########################################

##Use the normal function with the standard first RMassBankData-Input
##1), 2), 3) same deal, but the function for getting the spectra is different
#msmsRMBD <- newMsmsWorkspace()
#filesRMBD <- list.files(system.file("spectra", package="RMassBankData"),".mzML", full.names = TRUE)[1]
#msmsRMBD@files <- filesRMBD
#loadList(system.file("list/NarcoticsDataset.csv",package="RMassBankData"))
#msmsRMBD <- msmsWorkflow(msmsRMBD, mode="pH", steps=1)

##The Problem is that msmsXCMS@specs and msmsRMBD@specs should have the same format at the end
##Which means we need a function to convert the pseudospectrum from CAMERA to the list format that RMassBank uses
##I am attempting to create this function, although I have problems doing the specific parameters
##Mainly because there isn't nearly as much data in the pseudospectrum that findMsMsHRperxcms returns
##as there is in the specs that the original function returns

##msmsXCMS@specs is a list and has:
##"mz", "mzmin", "mzmax"
##"rt", "rtmin", "rtmax"
##"into", "intb", "maxo"     
##"sn", "isotopes", "adduct", "psg"
##Which are more or less just peaks with additional information

##msmsRMBD@specs is a list and has:
##$foundOK, which apparently checks if the data is ok? Something like that.
##$parentScan, which gives the parent a scan number? What is the parent supposed to be?
##$parentHeader, which is a vector containing 10 values for the parentscan
##$childScans, which gives the children their scan numbers
##$childHeader, which is a matrix containing 10 values for each of the childscans
##$parentPeak, which contains the peaks for the parent scan
##$peaks, which contains the peaks for the child scans
##$mz, which contains "mzMin", "mzMax", "mzCenter", where I don't know which maximum/minimum/center of what exactly it is supposed to describe.
##$id, which is the cpdid
##$formula, which is the formula
