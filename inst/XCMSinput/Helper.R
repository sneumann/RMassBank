##Helperscript to test a few things.
##This should work on any PC

library(RMassBank)
############################
############################
##CURRENTLY TRYING TO GET THIS TO WORK:
############################
############################

## xr <- system.file("microtofq/MSMSpos20_6.mzML", package="msdata", includeMsn=TRUE)
  
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
  
  
  
  
 handToRMB <- function(handSpecs, cpdID, MS1 = NA){
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
	ret$childHeader[1,4] <- length(handSpecs[,1])
	ret$childHeader[1,5] <- 0 ##Does this matter?
	ret$childHeader[1,6] <- findRt(cpdID)$RT * 60
	ret$childHeader[1,7] <- 0 ##Isn't possible per hand?
	ret$childHeader[1,8] <- 0 ##Isn't possible per hand?
	ret$childHeader[1,9] <- 0 ##Isn't possible per hand.
	ret$childHeader[1,10] <- 0 ##Isn't possible per hand.
	ret$childHeader[1,11] <- min(handSpecs[,1])
	ret$childHeader[1,12] <- max(handSpecs[,1])
	ret$childHeader[1,13] <- 1
	ret$childHeader[1,14] <- findMz(cpdID)[[3]]
	ret$childHeader[1,15] <- 1 ##Will be changed for different charges
	ret$childHeader[1,16] <- 0 ##There sadly isnt any precursor intensity to find in the msms-scans. Workaround?
	ret$childHeader[1,17:20] <- 0 ##Will be changed if merge is wanted
	ret$childHeader <- as.data.frame(ret$childHeader)
	ret$parentPeak <- matrix(nrow = 1, ncol = 2)
	colnames(ret$parentPeak) <- c("mz","int")
	ret$parentPeak[1,] <- c(findMz(cpdID)$mzCenter,100)
	ret$peaks <- list()
	ret$peaks[[1]] <- matrix(nrow = length(handSpecs[,1]), ncol = 2)
	colnames(ret$peaks[[1]]) <- c("mz","int")
	ret$peaks[[1]][,1] <- handSpecs[,1]
	ret$peaks[[1]][,2] <- handSpecs[,2]
	ret$mz <- findMz(cpdID)
	ret$id <- cpdID
	ret$formula <- findFormula(cpdID)
 }
  
 XCMStoRMB <- function(msmsXCMSspecs, cpdID, MS1 = NA){
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
	ret$parentPeak[1,] <- c(findMz(cpdID)$mzCenter,100)
	ret$peaks <- list()
	ret$peaks[[1]] <- matrix(nrow = length(msmsXCMSspecs[,1]), ncol = 2)
	colnames(ret$peaks[[1]]) <- c("mz","int")
	ret$peaks[[1]][,1] <- msmsXCMSspecs[,1]
	ret$peaks[[1]][,2] <- msmsXCMSspecs[,7]
	ret$mz <- findMz(cpdID)
	ret$id <- cpdID
	ret$formula <- findFormula(cpdID)
	return(ret)
}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
##This uses an already created template for the settings
##Change path or settings at own risk
loadRmbSettings(system.file("XCMSinput/mysettings.ini",package="RMassBank"))

##Use the xcms-CAMERA-Peakpicker with the ChelidonineMSn-mzData
##Extract filepath
##Load Compoundlist
##Get the spec into specsXCMS
msmsXCMS <- newMsmsWorkspace()
filesXCMS <- system.file("XCMSinput/Chelidonine_666_pos.mzData",package="RMassBank") 
msmsXCMS@files <- filesXCMS
loadList(system.file("XCMSinput/Chelidonine.csv",package="RMassBank"))

##Set the options correctly
rmbo <- getOption("RMassBank")
rmbo$annotations$entry_prefix <- 'IH'
rmbo$spectraList <- list(
  list(mode="CID",
       ces="20eV",
       ce="20eV",
       res=7500)
  )
#rmbo$recalibrateMS1 <- "none"
options("RMassBank" = rmbo)


#####WORKFLOW STEPS 1 to 8

########
##STEP 1
########

msmsXCMSspecs <- findMsMsHRperxcms.direct(msmsXCMS@files[1])
msmsXCMS@specs[[1]] <- XCMStoRMB(msmsXCMSspecs,666)
names(msmsXCMS@specs) <- basename(msmsXCMS@files)

########
##STEP 2 
########

mode = "pH"
msmsXCMS@analyzedSpecs <- lapply(msmsXCMS@specs, function(spec) {
				  s <- analyzeMsMs(spec, mode=mode, detail=TRUE, run="preliminary" )
				  return(s)
			  })
			  
########
##STEP 3
########

msmsXCMS@aggregatedSpecs <- aggregateSpectra(msmsXCMS@analyzedSpecs)

########
##STEP 4
########
recal <- makeRecalibration(msmsXCMS@aggregatedSpecs, mode)
msmsXCMS@rc <- recal$rc
msmsXCMS@rc.ms1 <- recal$rc.ms1
msmsXCMS@recalibratedSpecs <- recalibrateSpectra(mode, msmsXCMS@specs, w = msmsXCMS)

# p <- as.data.frame(msmsXCMS@specs[[1]]$peaks)
# rc <- msmsXCMS@rc
# Fix the column names so our
# prediction functions choose the right
# rows. 
# colnames(p) <- c("mzFound", "int")
# drecal <- predict(rc, newdata= p)
# Problem: Too many warnings

########
##STEP 5
########

msmsXCMS@analyzedRcSpecs <- lapply(msmsXCMS@recalibratedSpecs, function(spec){
		s <- analyzeMsMs(spec, mode=mode, detail=TRUE, run="recalibrated", cut=0, cut_ratio=0 )
		return(s)
	})
 for(f in msmsXCMS@files)
msmsXCMS@analyzedRcSpecs[[basename(as.character(f))]]$name <- basename(as.character(f))

########
##STEP 6
########

msmsXCMS@aggregatedRcSpecs <- aggregateSpectra(msmsXCMS@analyzedRcSpecs, addIncomplete=TRUE)
msmsXCMS@aggregatedRcSpecs$peaksUnmatchedC <- cleanElnoise(msmsXCMS@aggregatedRcSpecs$peaksUnmatched)

########
##STEP 7
########

msmsXCMS@reanalyzedRcSpecs <- reanalyzeFailpeaks(msmsXCMS@aggregatedRcSpecs, custom_additions="N2O", mode=mode)

########
##STEP 8
########

## This hack will not be neccessary after RMassBank 1.1.2
msmsXCMS@refilteredRcSpecs <- filterMultiplicity(msmsXCMS@reanalyzedRcSpecs, archivename=NA, mode)
msmsXCMS@refilteredRcSpecs$peaksOK <- msmsXCMS@refilteredRcSpecs$peaksFiltered
msmsXCMS@refilteredRcSpecs$peaksReanOK <- msmsXCMS@refilteredRcSpecs$peaksFilteredReanalysis

##
## Aim to be able to 
##  w <- msmsWorkflow(w, mode="pH", steps=c(2:8), archivename="pH_130301_pos")
##

############
##MBWORKFLOW
############

mbXCMS <- newMbWorkspace(msmsXCMS)
mbXCMS <- loadInfolists(mbXCMS, system.file("XCMSinput/infolists", package = "RMassBank"))
mbXCMS <- mbWorkflow(mbXCMS)
