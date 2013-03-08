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

##Use the xcms-CAMERA-Peakpicker with the Glucolesquerellin-MSn-mzData
##Extract filepath
##Load Compoundlist
##Get the spec into specsXCMS
msmsList <- newMsmsWorkspace()
fileList <- list.files(system.file("XCMSinput", package = "RMassBank"), "Glucolesquerellin", full.names=TRUE)[2:5]
msmsList@files <- fileList
loadList(system.file("XCMSinput/Chelidonine.csv",package="RMassBank"))

##Set the options right
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


#####WORKFLOW STEPS 1 to 8

########
##STEP 1
########

#msmsListspecs <- findMsMsHRperxcms.workflow(msmsList@files, mode="mH", method="centWave", peakwidth=c(5,10),
#												prefilter=c(3,200), ppm=25, snthr=5)
#msmsList@specs <- lapply(msmsList,XCMStoRMB, 2184)
#names(msmsList@specs) <- findName(2184)


#################
##STEP 1 MANUALLY, BECAUSE IT CAN'T FIND A FITTING PSEUDOSPECTRUM OF THE SECOND FILE?
#################
	fileName <- msmsList@files[2] 
	splitfn <- strsplit(fileName,'_')
    splitsfn <- splitfn[[1]]
    cpdID <- as.numeric(splitsfn[[length(splitsfn)-1]])
	
	parentMass <- findMz(cpdID)$mzCenter
	RT <- findRt(cpdID)$RT * 60
	mzabs <- 0.1
	
	getRT <- function(xa) {
		rt <- sapply(xa@pspectra, function(x) {median(peaks(xa@xcmsSet)[x, "rt"])})
	}
	##
	## MS
	##
	
	##
	## MSMS
	##
	xrmsms <- xcmsRaw(fileName, includeMSn=TRUE)

	## Where is the wanted isolation ?
	precursorrange <- range(which(xrmsms@msnPrecursorMz == parentMass)) ## TODO: add ppm one day

	## Fake MS1 from MSn scans
	xrmsmsAsMs <- msn2xcms(xrmsms)

	## Fake s simplistic xcmsSet
	xsmsms <-  xcmsSet (files=fileName,
					method="MS1")

	peaks(xsmsms) <- findPeaks(xrmsmsAsMs, method="centWave", peakwidth=c(5,10),
                         prefilter=c(3,200), ppm=25,
                         snthr=5, verbose.columns=T)

	## Get pspec 
	pl <- peaks(xsmsms)[,c("mz", "rt")]
	candidates <- which( pl[,"mz"] < parentMass + mzabs & pl[,"mz"] > parentMass - mzabs
						& pl[,"rt"] < RT * 1.1 & pl[,"rt"] > RT * 0.9 )
	
	anmsms <- xsAnnotate(xsmsms)
	anmsms <- groupFWHM(anmsms)

	## Now find the pspec for Chelidonine
	psp <- which(sapply(anmsms@pspectra, function(x) {candidates %in% x}))
	
	## Alternative: Spectrum closest to MS1
	##psp <- which.min(getRT(anmsms) - actualRT)