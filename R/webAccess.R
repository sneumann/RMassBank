#' @import XML RCurl rjson
NULL
## library(XML)
## library(RCurl)



#' Retrieve information from Cactus
#' 
#' Retrieves information from the Cactus Chemical Identifier Resolver
#' (PubChem).
#' 
#' It is not necessary to specify in which format the \code{identifier} is.
#' Somehow, cactus does this automatically.
#' 
#' @usage getCactus(identifier, representation)
#' @param identifier Any identifier interpreted by the resolver, e.g. an InChI
#' key or a SMILES code.
#' @param representation The desired representation, as required from the
#' resolver. e.g. \code{stdinchikey}, \code{chemspider_id}, \code{formula}...
#' Refer to the webpage for details.
#' @return The result of the query, in plain text. Can be NA, or one or
#' multiple lines (character array) of results.
#' @note Note that the InChI key is retrieved with a prefix (\code{InChIkey=}),
#' which must be removed for most database searches in other databases (e.g.
#' CTS).
#' @author Michael Stravs
#' @seealso \code{\link{getCtsRecord}}, \code{\link{getPcId}}
#' @references cactus Chemical Identifier Resolver:
#' \url{http://cactus.nci.nih.gov/chemical/structure}
#' @examples
#' 
#' # Benzene:
#' getCactus("C1=CC=CC=C1", "cas")
#' getCactus("C1=CC=CC=C1", "stdinchikey")
#' getCactus("C1=CC=CC=C1", "chemspider_id")
#' 
#' @export 
getCactus <- function(identifier, representation)
{

  ret <- tryCatch(
    getURLContent(paste(
      "http://cactus.nci.nih.gov/chemical/structure/",
      URLencode(identifier), "/", representation, sep='')),
    error = function(e) NA)
  if(is.na(ret))
    return(NA)
  if(ret=="<h1>Page not found (404)</h1>\n")
    return(NA)
  return(unlist(strsplit(ret, "\n")))
}

#' Search Pubchem CID
#' 
#' Retrieves PubChem CIDs for a search term.
#' 
#' Only the first result is returned currently. \bold{The function should be
#' regarded as experimental and has not thoroughly been tested.}
#' 
#' @usage getPcId(search)
#' @param search The search term.
#' @return The PubChem CID (in string type).
#' @author Michael Stravs
#' @seealso \code{\link{getCtsRecord}}, \code{\link{getCactus}}
#' @references PubChem search: \url{http://pubchem.ncbi.nlm.nih.gov/} 
#' 
#' Entrez E-utilities: \url{http://www.ncbi.nlm.nih.gov/books/NBK25500/}
#' @examples
#' 
#' # Benzene (again):
#' # Currently broken: getPcId("benzene")
#' 
#' @export
getPcId <- function(search)
{
  baseUrl <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
  term <- paste(baseUrl, "esearch.fcgi?db=pccompound&term=", URLencode(search), sep='')
  ret <-  getURL(term, timeout=5)
  errorvar <- 0
  currEnvir <- environment()
  tryCatch(
	{test <- getURL("www.google.com:81", timeout=5)},
	error=function(e){
	currEnvir$errorvar <- 1
	})
  
  if(errorvar){
	stop("Currently can't connect to pubchem")
  }
  
  #ret <- paste(ret, collapse='')
  s
  xml <- xmlParseDoc(ret,asText=TRUE)
  idNodes <- getNodeSet(xml, "/eSearchResult/IdList/Id")
  
  id <- xmlValue(idNodes[[1]])
  return(id)
}

# The following function is unfinished.
# getPcRecord <- function(pcid)
# {
#   baseUrl <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
#   term <- paste(baseUrl, "esummary.fcgi?db=pccompound&id=", URLencode(as.character(pcid)), 
#                 
#                 sep='')
#   ret <- getURL(term)
#   xml <- xmlParseDoc(ret,asText=TRUE)
#   browser()
# }


# Note: some of the CHEBI codes returned are erroneous. (When the entry in 
# CTS starts with "CHEBI:" instead of just the number, the XML output)
# Also, there is no ChemSpider ID in the XML output, unfortunately.








#' Retrieve information from CTS
#' 
#' Retrieves a complete CTS record from the InChI key.
#' 
#' @usage getCtsRecord(key)
#' 
#' @param key The InChI key. 
#' @return Returns a list with all information from CTS: \code{inchikey, 
#' 	inchicode, formula, exactmass} contain single values. \code{synonyms} contains
#' an unordered list of scored synonyms (\code{type, name, score}, where \code{type}
#' indicates either a normal name or a specific IUPAC name, see below).
#'  \code{externalIds} contains an unordered list of identifiers of the compound in 
#' various databases (\code{name, value}, where \code{name} is the database name and
#' \code{value} the identifier in that database.)
#' 
#' @note Currently, the CTS results are still incomplete; the name scores are all 0,
#' formula and exact mass return zero.
#' @references Chemical Translation Service:
#' \url{http://cts.fiehnlab.ucdavis.edu}
#' 
#' @examples
#' \dontrun{
#' data <- getCtsRecord("UHOVQNZJYSORNB-UHFFFAOYSA-N")
#' # show all synonym "types"
#' types <- unique(unlist(lapply(data$synonyms, function(i) i$type)))
#' }
#' \dontrun{print(types)}
#' 
#' @author Michele Stravs, Eawag <stravsmi@@eawag.ch>
#' @export
getCtsRecord <- function(key)
{
	baseURL <- "http://cts.fiehnlab.ucdavis.edu/service/compound/"
	data <- getURL(paste0(baseURL,key))
	r <- fromJSON(data)
	if(length(r) == 1)
		if(r == "You entered an invalid InChIKey")
			return(list())
	return(r)
}

#' Convert a single ID to another using CTS.
#' 
#' @usage getCtsKey(query, from = "Chemical Name", to = "InChIKey")
#' @param query ID to be converted
#' @param from Type of input ID
#' @param to Desired output ID 
#' @return An unordered array with the resulting converted key(s). 
#' 
#' @examples 
#' \dontrun{	
#' k <- getCtsKey("benzene", "Chemical Name", "InChIKey")
#' }
#' @author Michele Stravs, Eawag <stravsmi@@eawag.ch>
#' @export
getCtsKey <- function(query, from = "Chemical Name", to = "InChIKey")
{
	baseURL <- "http://cts.fiehnlab.ucdavis.edu/service/convert"
	url <- paste(baseURL, from, to, query, sep='/')
	data <- getURL(URLencode(url), timeout=7)
	r <- fromJSON(data)
	if(length(r) == 0)
		return(FALSE)
	else
	{
		# read out the results in simplest form:
		results <- unlist(lapply(r, function(row) row$result))
		return(results)
	}
}

#' Select a subset of external IDs from a CTS record.
#' 
#' @usage CTS.externalIdSubset(data, database)
#' @param data The complete CTS record as retrieved by \code{\link{getCtsRecord}}. 
#' @param database The database for which keys should be returned. 
#' @return Returns an array of all external identifiers stored in the record for the
#' given database.
#' 
#' @examples 
#' 
#' \dontrun{
#' # Return all CAS registry numbers stored for benzene.
#' data <- getCtsRecord("UHOVQNZJYSORNB-UHFFFAOYSA-N")
#' cas <- CTS.externalIdSubset(data, "CAS")
#' } 
#' 
#' @author Michele Stravs, Eawag <stravsmi@@eawag.ch>
#' @export
CTS.externalIdSubset <- function(data, database)
{
	select <- which(unlist(lapply(data$externalIds, function(id)
							{
								id[["name"]] == database
							})))
	keyEntries <- data$externalIds[select]
	keys <- unlist(lapply(keyEntries, function(e) e[["value"]]))
}

#' Find all available databases for a CTS record
#' 
#' @usage CTS.externalIdTypes(data)
#' @param data The complete CTS record as retrieved by \code{\link{getCtsRecord}}.  
#' @return Returns an array of all database names for which there are external 
#' identifiers stored in the record.
#' 
#' @examples 
#' 
#' \dontrun{
#' # Return all databases for which the benzene entry has
#' # links in the CTS record.
#' 
#' data <- getCTS("UHOVQNZJYSORNB-UHFFFAOYSA-N")
#' databases <- CTS.externalIdTypes(data)
#' } 
#' 
#' @author Michele Stravs, Eawag <stravsmi@@eawag.ch>
#' @export
CTS.externalIdTypes <- function(data)
{
	unique(unlist(lapply(data$externalIds, function(id)
							{
								id[["name"]]
							})))
}
