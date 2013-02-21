##Helperscript to test a few things.
##This should work on any PC

library(RMassBank)

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
msmsXCMS@specs <- findMsMsHRperxcms(msmsXCMS@files[1]) ##

##Use the normal function with the standard first RMassBankData-Input
##1), 2), 3) same deal, but the function for getting the spectra is different
msmsRMBD <- newMsmsWorkspace()
filesRMBD <- list.files(system.file("spectra", package="RMassBankData"),".mzML", full.names = TRUE)[1]
msmsRMBD@files <- filesRMBD
loadList(system.file("list/NarcoticsDataset.csv",package="RMassBankData"))
msmsRMBD <- msmsWorkflow(msmsRMBD, mode="pH", steps=1)

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

##The problematic parts are the headers,mzMin/-Max/-Center and the differentiation between childscans and parentscans.
##I can just look up the rest in the compoundlist