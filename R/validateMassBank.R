# The function that validates a MassBank record
# 
#' MassBank record validator
#' 
#' Validates MassBank records in plaintext-format 
#' 
#' 
#'
#' @param File A path to the records that should be validated.
#' @return Nothing
#' @seealso \code{\link{parseMassBank}}
#' @author Erik Mueller
#' @examples \dontrun{
#' 		validate(c("filepath_to_records/RC00001.txt","filepath_to_records/RC00002.txt))
#' }
#' @export
validate <- function(File){
	# Is the argument a directory?
	# If yes, list the files
	RMassBank.env$Instrument_List <- .getInstruments()
	if(file.info(File[1])$isdir){
	    Files <- list.files(path = File, full.names = TRUE)
	}

	# Parsing with the help the parseMassBank-function
	RMassBank.env$mb <- lapply(Files,parseMassBank)
	print("yeah")
	# Test RMassBank Objects with RUnit
	# This loop creates the tests and defines one test suite for every record
	tests <- list()
	for(i in 1:length(Files)){
		if(RMassBank.env$mb[[i]]@compiled_ok[[1]][['AC$MASS_SPECTROMETRY']][['MS_TYPE']] == "MS2" || RMassBank.env$mb[[i]]@compiled_ok[[1]][['AC$MASS_SPECTROMETRY']][['MS_TYPE']] == "MS"){
		tests[[i]] <- defineTestSuite(Files[i], dirs = system.file(package="RMassBank", "unitTests"), testFileRegexp = "^runit.MS2.test.[rR]$",
                testFuncRegexp = "^test.+",
                rngKind = "Marsaglia-Multicarry",
                rngNormalKind = "Kinderman-Ramage")
		} else{
			tests[[i]] <- defineTestSuite(Files[i], dirs = system.file(package="RMassBank", "unitTests"), testFileRegexp = "^runit.MSn.test.[rR]$",
                testFuncRegexp = "^test.+",
                rngKind = "Marsaglia-Multicarry",
                rngNormalKind = "Kinderman-Ramage")
		}
	}
	print("Starting Tests")
	# Testing the list of Testsuites
	testData <- runTestSuite(tests)
	# Prints the HTML-record
	printHTMLProtocol(testData, fileName = paste(system.file(package = "RMassBank"),"report.html", sep = ""))
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

# This function uses the SMILES-formula to calculate the mass
# It is needed for a test and needs to be exported so RUnit can use it
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