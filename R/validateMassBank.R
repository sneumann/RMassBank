#' Validate a record against a set of RUnit tests
#' 
#' Validates a plaintext-MassBankrecord, or recursively all records
#' below a directory. The Unit Tests to be used
#' are installed in RMassBank/
#' 
#' @aliases validate
#' @usage validate(path)
#' @param path The filepath to a single record, or a directory to search recursively
#' @examples
#' \dontrun{
#' validate("/tmp/MassBank/OpenData/record/")
#' }
#' @export
validate <- function(path) {

        if (!require(ontoCAT)) {
          error("Package ontoCAT missing. Validation requires package ontoCAT and RUnit")
        }

        if (!require(RUnit)) {
          error("Package RUnit missing. Validation requires package ontoCAT and RUnit")
        }

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
		tests[[i]] <- defineTestSuite(Files[i], dirs = system.file(package="RMassBank", "unitTests"), testFileRegexp = "runit.MS2.test.R",
                #testFuncRegexp = "^test.+",
                rngKind = "Marsaglia-Multicarry",
                rngNormalKind = "Kinderman-Ramage")
		} else{
			tests[[i]] <- defineTestSuite(Files[i], dirs = system.file(package="RMassBank", "unitTests"), testFileRegexp = "^runit.MSn.test.[rR]$",
                #testFuncRegexp = "^test.+",
                rngKind = "Marsaglia-Multicarry",
                rngNormalKind = "Kinderman-Ramage")
		}
	}
	print("Starting Tests")
	# Testing the list of Testsuites
	testData <- runTestSuite(tests)
	# Prints the HTML-record
	printHTMLProtocol(testData, fileName = paste(getwd(),"/report.html", sep = ""))
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
