##Helperscript to test a few things.
##This should work on any PC

library(RMassBank)
############################
############################
##CURRENTLY TRYING TO GET THIS TO WORK:
############################
############################

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
## getOption("RMassBank")$spectraList

## spectraList:
##  # First scan: CID 20
## - mode: CID
##   ces: 20
##   ce: 20
##   res: 7500

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
##Extract filepath
##Load Compoundlist
##Get the spec into specsXCMS
msmsXCMS <- newMsmsWorkspace()
filesXCMS <- system.file("XCMSinput/Chelidonine_666_pos.mzData",package="RMassBank") 
msmsXCMS@files <- filesXCMS
loadList(system.file("XCMSinput/Chelidonine.csv",package="RMassBank"))

#####WORKFLOW STEPS 1 to 8

########
##STEP 1
########

msmsXCMSspecs <- findMsMsHRperxcms.direct(msmsXCMS@files[1])
msmsXCMS@specs[[1]] <- XCMStoRMB(msmsXCMSspecs,666)
names(msmsXCMS@specs) <- findName(666)

########
##STEP 2 
########

msmsXCMS@specs[[1]]$parentPeak[,"mz"] <- 354.134499255997
msmsXCMS@specs[[1]]$parentPeak[,"int"] <- 354.134499255997

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
#recal <- makeRecalibration(msmsXCMS@aggregatedSpecs, mode)
#msmsXCMS$rc <- recal$rc
#msmsXCMS$rc.ms1 <- recal$rc.ms1
#msmsXCMS@recalibratedSpecs <- recalibrateSpectra(mode, msmsXCMS@specs, w = msmsXCMS)

recalibrate.addMS1data(spec = msmsXCMS@aggregatedSpecs, mode = "pH", 15)

rmbo <- getOption("RMassBank")
rmbo$spectraList <- list(
  list(mode="CID",
       ces="10eV",
       ce="10eV",
       res=12000),
  list(mode="CID",
       ces="20eV",
       ce="20eV",
       res=12000),
  list(mode="CID",
       ces="30eV",
       ce="30eV",
       res=12000),
  list(mode="CID",
       ces="40eV",
       ce="40eV",
       res=12000)
  )

options("RMassBank" = rmbo )

