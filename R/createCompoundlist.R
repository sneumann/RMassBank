#' @import readJDX
#' @import webchem
#' @import data.table
#' @import ChemmineR
#' @import ChemmineOB

#' @title Add a header to a Multiblock JCAMP file
#'
#' @description JCAMP files containing multiple blocks are usually structured
#' by so-called link blocks. If no link block is present, the readJDX
#' package is not able to parse the file. This method will add a link
#' block at the top of the given file or print a message if an existing
#' link block is found. The file is not changed in this case.
#'
#' @param filename character
#' The name of the file to which a link block should be added.
#' The filename is also used as content for the TITLE field in the link block
#' @return Nothing is returned
#' @examples \dontrun{
#' 			updateHeader("my_multiblock_jcamp.jdx")
#' }
#' @author pstahlhofen
#' @export
updateHeader <- function(filename) {
	lines <- readLines(filename)
	block_pattern <- "##BLOCKS=(.*)"
	contains_header <- any(grepl(block_pattern, lines))
	if (contains_header) {
		cat('Header is already present. No update performed\n')
	}
	else {
		end_pattern <- '##END='
		n_blocks <- sum(grepl(end_pattern, lines))
		field_names <- paste0('##', c('TITLE', 'BLOCKS', 'DATA TYPE'))
		field_values <- c(filename, n_blocks, 'LINK')
		header_block <- paste(field_names, field_values, sep='=')
		updated <- c(header_block, lines)
		writeLines(updated, filename)
		cat('Header block added successfully\n')
	}
}

#' Get the content of a field in a JCAMP file
#'
#' The content will always be returned as character-string
#'
#' @param parsedJDX list as created by readJDX
#' A parsed, single-block JCAMP file
#' @param field_name character
#' The name of the field (e.g. 'CAS REGISTRY NO')
#' @return The field's content
#' @examples \dontrun{
#' 			parsedJDX <- readJDX('my_singleblock_jcamp.dx')
#' 			title <- getField(parsedJDX, "TITLE")
#' }
#' @author pstahlhofen
#' @seealso readJDX
#' @export
getField <- function(parsedJDX, field_name) {
	field <- grep(field_name, parsedJDX$metadata, value=TRUE)
	field_split <- strsplit(field, '=')[[1]]
	field_value <- field_split[-1]
	return(field_value)
}

getCAS <- function(parsedJDX) {return(getField(parsedJDX, 'CAS REGISTRY NO'))}

getTitle <- function(parsedJDX) {return(getField(parsedJDX, 'TITLE'))}

#' Convert CAS to SMILES
#' 
#' This is a wrapper for \code{webchem::cir_query}, using the
#' CACTUS API at https://cactus.nci.nih.gov/chemical/structure_documentation
#' for the conversion. Before converting the CAS number, the 
#' name is checked whether it contains the word 'derivative'.
#' If so, the conversion is stopped and NA is returned.
#' Also, a warning will be printed in this case.
#'
#' The API allows only one query per second. This is a hard-
#' coded feature
#'
#' @param CAS_number character
#' The CAS registry number of a compound
#' @param name character
#' The compound's name
#' @return The SMILES code of the compound as character-string
#' @examples SMILES_ethanol <- CAS2SMILES("64-17-5", "Ethanol")
#' @author pstahlhofen
#' @export
CAS2SMILES <- function(CAS_number, name) {
	if(grepl('derivative', name)) {
		warning(paste("Converting CAS to SMILES for the compound",
			      name, "might yield a wrong result.",
			      "Please provide the structure manually.",
			      sep=" "))
		return(NA)
	}
	return(cir_query(CAS_number, from='cas', to='smiles'))
}

#' Create a Compoundlist from JCAMP files
#'
#' This method will automatically look for all single-block
#' JCAMP files in the directory by picking all files ending in '.dx'
#' (and not '.jdx'). A csv-file named 'Compoundlist.csv' will
#' be created in the same directory. The Compoundlist contains
#' columns 'ID', 'Name', 'SMILES' and 'CAS' where 'SMILES' might
#' be empty if the compound is a derivative or if the CAS number
#' could not be converted (see CAS2SMILES).
#'
#' @return This method has no return value.
#' @examples \dontrun{
#' 			# Prepare the compoundlist-creation
#' 			splitMultiblockDX('my_multiblock_jcamp.jdx')
#' 			createCompoundlist()
#' }
#' @author pstahlhofen
#' @seealso CAS2SMILES
#' @export
createCompoundlist <- function() {
	files <- list.files(getwd(), pattern='[^j]dx$')
	parsedFiles <- lapply(files, readJDX)
	CAS_numbers <- sapply(parsedFiles, getCAS)
	names <- sapply(parsedFiles, getTitle)
	SMILES_codes <- sapply(seq_along(names), function(idx) {
		return(CAS2SMILES(CAS_numbers[idx], names[idx]))
	})
	compoundlist <- data.frame(ID=seq_along(names),
				   Name=names,
				   SMILES=unlist(SMILES_codes),
				   CAS=CAS_numbers)
	fwrite(compoundlist, file='Compoundlist.csv')
}

#' Filter a Compoundlist for missing SMILES values
#'
#' Read the Compoundlist given by the filename and write a
#' 'Compoundlist_filtered.csv', containing only the lines
#' with a SMILES string
#'
#' @param filename character
#' The name of the csv-file to be read
#' @examples \dontrun{
#' 			filterCompoundlist('Compoundlist.csv')
#' }
#' @return This method has no return value.
#' @author pstahlhofen
#' @export
filterCompoundlist <- function(filename) {
	compoundlist <- fread(filename)
	filtered <- compoundlist[which(compoundlist$SMILES!=""), ]
	fwrite(filtered, file='Compoundlist_filtered.csv')
}

#' Convert a Compoundlist into an SDF
#'
#' The resulting SDF will be written to a file named 'Compoundlist.sdf'.
#' The header for each block is the chemical name, tags for ID, SMILES and CAS
#' are added in the description block
#'
#' @param filename character
#' The name of the csv-file to be read. Note that the compoundlist
#' has to be filtered already.
#' @return This method has no return value.
#' @examples \dontrun{
#' 			compoundlist2SDF("Compoundlist_filtered.csv")
#' }
#' @author pstahlhofen
#' @export
compoundlist2SDF <- function(filename) {
	compoundlist <- fread(filename)
	SMILES <- compoundlist$SMILES
	if (any(SMILES=="")) {
		stop(paste("The provided compoundlist must be filtered",
		  "for missing SMILES values first.", sep=" "))
	}
	names(SMILES) <- compoundlist$Name
	SDFset <- smiles2sdf(SMILES)
	valid <- validSDF(SDFset)
	if (!all(valid)) {
		invalid <- names(SMILES[!valid])
		warning_message <- paste('The following compounds',
		  'cannot be converted to SDF blocks:')
		warning(paste(c(warning_message, invalid), sep='\n\t- '))
	}
	SDFset <- SDFset[valid]
	SMILES <- SMILES[valid]
	IDs <- compoundlist[valid, ID]
	CAS <- compoundlist[valid, CAS]
	SDFset@SDF <- lapply(seq_along(SDFset), function(idx) {
		single_SDF <- SDFset[[idx]]
		metadata <- c(IDs[idx], SMILES[idx], CAS[idx])
		names(metadata) <- c('ID', 'SMILES', 'CAS')
		single_SDF@datablock <- metadata
		return(single_SDF)
	})
	write.SDF(SDFset, 'Compoundlist.sdf', cid=TRUE)
}

