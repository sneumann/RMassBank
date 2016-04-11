w <- newMsmsWorkspace()
RmbDefaultSettings()
files <- list.files(system.file("spectra", package="RMassBankData"), ".mzML", full.names = TRUE)

w@files <- files
loadList(system.file("list/NarcoticsDataset.csv", package="RMassBankData"))
w <- msmsWorkflow(w, mode="pH", steps=c(1:4), archivename="pH_narcotics")
w <- msmsWorkflow(w, mode="pH", steps=c(5:8), archivename="pH_narcotics")	
mb <- newMbWorkspace(w)
mb <- loadInfolists(mb, system.file("infolists", package="RMassBankData"))
mb <- mbWorkflow(mb)


### Error message code structure borrowed from rcppgls:
### https://github.com/eddelbuettel/rcppgsl/blob/master/inst/unitTests/runTests.R

testElec <- RUnit::defineTestSuite("Electronic noise and formula calculation Test", dirs = system.file("unitTests", 
									package="RMassBank"), testFileRegexp = "test_El.*")
testAcqu <- RUnit::defineTestSuite("Evaluation of data acquisition process", dirs = system.file("unitTests", 
									package="RMassBank"), testFileRegexp = "test_mzR.R")
testNope <- RUnit::defineTestSuite("Evaluation of correct handling if no peaks are found", dirs = system.file("unitTests", 
									package="RMassBank"), testFileRegexp = "test_NOPEAKS.R")

errEval <- function(testSuite){
	Res <- RUnit::runTestSuite(testSuite)
	err <- getErrors(Res)
	if((err$nFail + err$nErr) > 0){
            stop(sprintf('unit test problems: %d failures, %d errors in testfunction: "%s"', err$nFail, err$nErr, testSuite$name) )
        } else{
            success <- err$nTestFunc - err$nFail - err$nErr - err$nDeactivated
            cat( sprintf( "%d / %d\n", success, err$nTestFunc ) )
    }
}

errEval(testElec)
errEval(testAcqu)
errEval(testNope)