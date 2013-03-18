##Helperscript to test a few things.
##This should work on any PC

library(RMassBank)
  
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 
  ################################################################################## 


##This uses an already created template for the settings
##Change path or settings at own risk
loadRmbSettings(system.file("XCMSinput/mysettings.ini",package="RMassBank"))

##Use the xcms-CAMERA-Peakpicker with the Glucolesquerellin-MSn-mzData
##Extract filepath
##Load Compoundlist
##Get the spec into specsXCMS
msmsList <- newMsmsWorkspace()
fileList <- list.files(system.file("XCMSinput", package = "RMassBank"), "Glucolesquerellin", full.names=TRUE)[2:3]
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
       res=12000)#,
#  list(mode="CID",
#       ces="40eV",
#       ce="40eV",
#       res=12000)
  )
options("RMassBank" = rmbo)

#####WORKFLOW STEPS 1 to 8

########
##STEP 1
########

print("Step 1")
handSpecs <- matrix(0,4,2)
handSpecs[,1] <- c(274.986685367956, 259.012401087427, 95.9493025990907, 96.9573002472772)
handSpecs[,2] <- c(357,761, 2821, 3446)
hand <- list()
hand[[1]] <- handSpecs
mode="mH"
msmsList <- findMsMsHRperX.workflow(msmsList, mode="mH", method="centWave", peakwidth=c(5,10),
												prefilter=c(3,200), ppm=25, snthr=5)
#msmsList <- addHand(msmsList,hand)
#msmsList <- addMB(msmsList,"record/Fukuyama_Univ/FU000001.txt")

########
##STEP 2 
########
print("Step 2")
mode = "mH"
msmsList@analyzedSpecs <- lapply(msmsList@specs, function(spec) {
				  s <- analyzeMsMs(spec, mode=mode, detail=TRUE, run="preliminary" )
				  return(s)
			  })
			  
########
##STEP 3
########
print("Step 3")
msmsList@aggregatedSpecs <- aggregateSpectra(msmsList@analyzedSpecs)

########
##STEP 4
########
print("Step 4")
recal <- makeRecalibration(msmsList@aggregatedSpecs, mode)
msmsList@rc <- recal$rc
msmsList@rc.ms1 <- recal$rc.ms1
msmsList@recalibratedSpecs <- recalibrateSpectra(mode, msmsList@specs, w = msmsList)

# p <- as.data.frame(msmsList@specs[[1]]$peaks)
# rc <- msmsList@rc
# Fix the column names so our
# prediction functions choose the right
# rows. 
# colnames(p) <- c("mzFound", "int")
# drecal <- predict(rc, newdata= p)
# Problem: Too many warnings

########
##STEP 5
########
print("Step 5")
msmsList@analyzedRcSpecs <- lapply(msmsList@recalibratedSpecs, function(spec){
		s <- analyzeMsMs(spec, mode=mode, detail=TRUE, run="recalibrated", cut=0, cut_ratio=0 )
		return(s)
	})
 for(f in msmsList@files)
msmsList@analyzedRcSpecs[[basename(as.character(f))]]$name <- basename(as.character(f))

########
##STEP 6
########
print("Step 6")
msmsList@aggregatedRcSpecs <- aggregateSpectra(msmsList@analyzedRcSpecs, addIncomplete=TRUE)
msmsList@aggregatedRcSpecs$peaksUnmatchedC <- cleanElnoise(msmsList@aggregatedRcSpecs$peaksUnmatched)

########
##STEP 7
########
print("Step 7")
msmsList@reanalyzedRcSpecs <- reanalyzeFailpeaks(msmsList@aggregatedRcSpecs, custom_additions="N2O", mode=mode)

########
##STEP 8
########
print("Step 8")
msmsList@refilteredRcSpecs <- filterMultiplicity(msmsList@reanalyzedRcSpecs, archivename=NA, mode)
msmsList@refilteredRcSpecs$peaksOK <- msmsList@refilteredRcSpecs$peaksFiltered
msmsList@refilteredRcSpecs$peaksReanOK <- msmsList@refilteredRcSpecs$peaksFilteredReanalysis

############
##MBWORKFLOW
############

mbList <- newMbWorkspace(msmsList)
mbList <- loadInfolists(mbList, system.file("XCMSinput/infolists2", package = "RMassBank"))
mbList <- mbWorkflow(mbList)
