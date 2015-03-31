#' Validate MassBank records with a set of Unit tests
#' 
#' Validates a plain text MassBank record, or recursively all
#' records within a directory. The Unit Tests to be used are
#' installed in RMassBank/inst/validationTests and currently include 
#' checks for NAs, peaks versus precursor, precursor mz, 
#' precursor type, SMILES vs exact mass, total intensities and
#' title versus type. The validation report is saved as 
#' "report.html" in the working directory.
#' 
#' @aliases validate
#' @usage validate(path, simple = TRUE)
#' @param path The filepath to a single record, or a directory to search recursively
#' @param simple If TRUE the function creates a simpler form of the RUnit .html report, better readable for humans. If FALSE it returns the unchanged RUnit report.
#' @examples
#' \dontrun{
#' validate("/tmp/MassBank/OpenData/record/")
#' }
#' @export
validate <- function(path, simple = TRUE) {

        require(ontoCAT)
		require(RUnit)

	# Is the argument a directory?
	# If yes, list the files
	RMassBank.env$Instrument_List <- .getInstruments()
	RMassBank.env$testnumber <- 1
	if(file.info(path[1])$isdir){
	    Files <- list.files(path = path,
                                recursive=TRUE, 
                                full.names = TRUE)
	} else {Files <- path}
	# Parsing with the help the parseMassBank-function
	RMassBank.env$mb <- lapply(Files,parseMassBank)
	# Test RMassBank Objects with RUnit
	# This loop creates the tests and defines one test suite for every record
	tests <- list()
	for(i in 1:length(RMassBank.env$mb)){
		if(RMassBank.env$mb[[i]]@compiled_ok[[1]][['AC$MASS_SPECTROMETRY']][['MS_TYPE']] == "MS2" || RMassBank.env$mb[[i]]@compiled_ok[[1]][['AC$MASS_SPECTROMETRY']][['MS_TYPE']] == "MS"){
		tests[[i]] <- defineTestSuite(Files[i], dirs = system.file(package="RMassBank", "validationTests"), testFileRegexp = "runit.MS2.test.R",
                #testFuncRegexp = "^test.+",
                rngKind = "Marsaglia-Multicarry",
                rngNormalKind = "Kinderman-Ramage")
		} else{
			tests[[i]] <- defineTestSuite(Files[i], dirs = system.file(package="RMassBank", "validationTests"), testFileRegexp = "^runit.MSn.test.[rR]$",
                #testFuncRegexp = "^test.+",
                rngKind = "Marsaglia-Multicarry",
                rngNormalKind = "Kinderman-Ramage")
		}
	}
	print("Starting Tests")
	# Testing the list of Testsuites
	testData <- suppressWarnings(runTestSuite(tests,verbose=0))
	# Prints the HTML-record
	printHTMLProtocol(testData, fileName = paste0(getwd(),"/report.html"))
	if(simple){
		fileConnection <- file(paste0(getwd(),"/report.html"), open = "r")
		htmlFile <- readLines(fileConnection)
		close(fileConnection)
		htmlFile <- gsub(">test.NA", ">No NAs contained in peak list", htmlFile)
		htmlFile <- gsub(">test.peaksvsprecursor", ">One peak m/z  with no noticable difference from the precursor mass", htmlFile)
		htmlFile <- gsub(">test.precursormz", ">Mass of precursor m/z possible with given mass and type", htmlFile)
		htmlFile <- gsub(">test.PrecursorType", ">Precursor type valid", htmlFile)
		htmlFile <- gsub(">test.smilesvsexactmass", ">Smiles code represents a molecule with specified exact mass", htmlFile)
		htmlFile <- gsub(">test.sumintensities", ">All intensies greater than zero", htmlFile)
		htmlFile <- gsub(">test.TitleVsType", ">Precursor type are the same in the title and in the document", htmlFile)
		htmlFile <- gsub("\\(1 checks\\)", "", htmlFile)
		htmlFile <- gsub("\\({1}.{1,6}seconds\\)", "", htmlFile)
		htmlFile <- gsub("Test Suite: ", "", htmlFile)
		htmlFile <- gsub("h5", "h2", htmlFile)
		##Remove ending
		poshr <- grep("<hr>", htmlFile, fixed=TRUE)
		poshr <- poshr[length(poshr)]
		htmlFile <- htmlFile[1:(poshr)]
		
		ullines <- grep("<ul>", htmlFile)
		htmlFile[ullines] <- gsub("</a><ul>", "</a></li>", htmlFile[ullines])
		htmlFile[ullines] <- gsub("</li></ul></li></ul>", "</li></ul>", htmlFile[ullines])
		htmlFile[ullines] <- gsub('<li><a href=".+">Test file: runit.MS2.test.R</a></li>', "", htmlFile[ullines])
		htmlFile[ullines] <- gsub('</a>.+RMassBank/validationTests<br/>',"</a>", htmlFile[ullines])
		
		##Remove superfluous information
		htmlFile <- gsub("Test function regexp: ^test.+<br/>Test file regexp: runit.MS2.test.R<br/>Involved directory:<br/>", "", htmlFile, fixed=TRUE)
		
		fileConnection <- file(paste0(getwd(),"/report.html"), open = "w")
		writeLines(htmlFile, fileConnection)
		close(fileConnection)
	}
	print(paste("Report for the file(s) finished"))
}

# This function checks if an .obo-file is readable for ontoCAT
.isOboReadable <- function(filename){

	# getOntology() has a problem with reading relative Windows paths(it wants an URI),
	# so the path has to be made absolute
	# I reckon this should work under Linux without doing that
	ont <- getOntology(normalizePath(filename))
	if(is.null(getOntologyAccession(ont))){
		return(FALSE)
	}
	return(TRUE)
}

# This function downloads the psi-ms.obo-ontology so we can get the allowed instrument-names
# This is a _temporary_ fix until I find out why getOntology() doesn't work when there are "import:"-lines in the .obo-file
# Until then I will simply remove them, because we don't need the imported ontologies
.downloadPsiObo <- function(){
		connPsiObo <- url("http://psidev.cvs.sourceforge.net/viewvc/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo")
		oboFile <- readLines(connPsiObo)
		close(connPsiObo)
		oboFile <- oboFile[-grep("import:",oboFile)] 
		connLocal <- file("psi-ms.obo")
		writeLines(oboFile,connLocal)
		close(connLocal)
}

# Checks if the psi-ms.obo is there
# Will be converted to "checkforinstruments" as soon as I can find the problem
# with getOntology()
.checkForPsiMs <- function(){
	
	if(file.exists("psi-ms.obo")){
		if(.isOboReadable("psi-ms.obo")){
			print("It seems that you have a working psi-ms.obo, do you want to update it? [y/n]")
			while(TRUE){
				answer <- readLines(stdin(), n=1, warn=FALSE)
				if(answer == "y"){
					.downloadPsiObo()
					return(TRUE)
				}
				if(answer == "n"){
					return(TRUE)
				}
				print("Please type exactly y or n")
			}
		}
	}
	.downloadPsiObo()
	return(TRUE)
}

# This is a list of the possible instrument names 
.getInstruments <- function(){
	Onto <- getOntology(system.file(package = "RMassBank", "psi-ms.obo"))
	instrumentTerms <- getAllTermChildrenById(Onto,"MS_1000031")
	instruments <- vector()	
	for(i in 1:length(instrumentTerms)){
		instruments[i] <- getLabel(instrumentTerms[[i]])
	}
	return(instruments)
}

#' Calculate the mass from a SMILES-String
#' 
#' Uses a SMILES-String to calculate the mass using rcdk-integrated functions.
#'
#' @aliases smiles2mass
#' @usage smiles2mass(SMILES)
#' @param SMILES A String-object representing a SMILES
#' @return The calculated mass of the given SMILES-Formula
#' @author Erik Mueller
#' @examples \dontrun{
#' 		smiles2mass("CC(=O)NC(C(O)1)C(O)C(OC(O2)C(O)C(OC(O3)C(O)C(O)C(O)C(CO)3)C(O)C(CO)2)C(CO)O1")
#' }
#' @export
smiles2mass <- function(SMILES){
	massfromformula <- parse.smiles(SMILES)[[1]]
	do.typing(massfromformula)
	do.aromaticity(massfromformula)
	convert.implicit.to.explicit(massfromformula)
	do.isotopes(massfromformula)
	mass <- get.exact.mass(massfromformula)
	return(mass)
}

.unitTestRMB <- function(WD=getwd()){
	require(RUnit)
	library(RMassBank)
	library(RMassBankData)
	oldwd <- getwd()
	setwd(WD)
	w <- newMsmsWorkspace()
	RmbDefaultSettings()
	files <- list.files(system.file("spectra", package="RMassBankData"),
		 ".mzML", full.names = TRUE)
	basename(files)
	# To make the workflow faster here, we use only 2 compounds:
	w@files <- files
	loadList(system.file("list/NarcoticsDataset.csv", 
		package="RMassBankData"))
	w <- msmsWorkflow(w, mode="pH", steps=c(1:4), archivename = 
					"pH_narcotics")
	w <- msmsWorkflow(w, mode="pH", steps=c(5:8), archivename = 
			"pH_narcotics")	
	w2 <- newMbWorkspace(w)
	#w2 <- mbWorkflow(w2)
	#w2 <- loadInfolist(w2, "infolist.csv")
	#w2 <- mbWorkflow(w2)
	
	testSuite <- defineTestSuite("Electronic noise and formula calculation Test", dirs = system.file("unitTests", 
		package="RMassBank"), testFileRegexp = "runit.EN_FC.R",
					#testFuncRegexp = "^test.+",
					rngKind = "Marsaglia-Multicarry",
					rngNormalKind = "Kinderman-Ramage")
					
	testData <- suppressWarnings(runTestSuite(testSuite))
	
	file.remove(c("pH_narcotics_Failpeaks.csv","pH_narcotics.RData","pH_narcotics_RA.RData","pH_narcotics_RF.RData"))
	
	# Prints the HTML-record
	printTextProtocol(testData)
	setwd(oldwd)
	return(testData)
}
