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
loadList(system.file("XCMSinput/compoundList.csv",package="RMassBank"))

msmsHand <- newMsmsWorkspace()
filesHand <- list.files(system.file("XCMSinput", package = "RMassBank"), "Glucolesquerellin.+(csv)$", full.names=TRUE)
msmsHand@files <- filesHand

msmsMB <- newMsmsWorkspace()
filesMB <- "C:/Users/Erik/Downloads/JP000001.txt"
msmsMB@files <- filesMB

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

msmsXCMS <- msmsWorkflow(msmsXCMS,steps=1, method="xcms")
msmsHand <- msmsWorkflow(msmsHand,steps=1, method="peaklist")
#msmsMB	 <- msmsWorkflow(msmsMB,steps=1,method="MassBank")

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
print("Step 3")
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
print("Step 5")
########
##STEP 6
########

msmsXCMS@aggregatedRcSpecs <- aggregateSpectra(msmsXCMS@analyzedRcSpecs, addIncomplete=TRUE)
msmsXCMS@aggregatedRcSpecs$peaksUnmatchedC <- cleanElnoise(msmsXCMS@aggregatedRcSpecs$peaksUnmatched)
print("Step 6")
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
