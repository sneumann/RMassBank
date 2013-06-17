#' @import XML RCurl
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
#' getPcId("benzene")
#' 
#' @export
getPcId <- function(search)
{
  
  baseUrl <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
  term <- paste(baseUrl, "esearch.fcgi?db=pccompound&term=", URLencode(search), sep='')
  ret <-  getURL(term)
  #ret <- paste(ret, collapse='')
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
#' Retrieves chemical information about a compound from Chemical Translation
#' Service (CTS) from a known identifier.
#' 
#' 
#' @usage getCtsRecord(key, from = "inchikey", to = c("cas", "hmdb", "kegg",
#' "sid", "chebi", "inchi", "lipidmap", "smiles", "cid", "inchikey", "mass",
#' "formula", "iupac", "names"))
#' @param key The search term (or key).
#' @param from The format of the key. Allowed are \code{"cas", "hmdb", "kegg",
#' "sid", "chebi", "inchi", "lipidmap", "smiles", "cid", "inchikey", "mass",
#' "formula", "iupac", "name"}.
#' @param to The list of result types which should be returned. Allowed are
#' \code{"cas", "hmdb", "kegg", "sid", "chebi", "inchi", "lipidmap", "smiles",
#' "cid", "inchikey", "mass", "formula", "iupac", "name"}.
#' @return Returns a named list with the values of the results. The list item
#' \code{"names"} is a matrix with columns \code{"name", "score"}, with
#' \code{score} being an indicator of the reliability of the name assignment.
#' @note The return values are not 100% reliable, e.g. a known bug returns
#' "ChEBI" for the \code{chebi} entry instead of the actual ChEBI code in some
#' instances.
#' @author Michael Stravs
#' @seealso \code{\link{getCactus}},\code{\link{getPcId}}
#' @references Chemical Translation Service:
#' \url{http://uranus.fiehnlab.ucdavis.edu:8080/cts/homePage}
#' @examples
#' 
#' getCtsRecord("benzene", "name")
#' 
#' @export
getCtsRecord <- function(key, from = "inchikey", 
  to = c("cas","hmdb","kegg","sid","chebi","inchi","lipidmap","smiles","cid",
         "inchikey","mass","formula","iupac","names"))
{
	require(devtools)
	require(RJSONIO)
	#install_github(repo = "CTSgetR", username = "dgrapov")
	require(CTSgetR)
	# checks
	if(from %in% c("", "None", "Unknown", "Not available"))
    return(NA)
	
	ChemSpID<-CTSgetR(key,from="InChIKey",to="ChemSpider",parallel=FALSE)
	CTS.options<-CTS.options()[2:10]
	CTS.options # see options
	id<-ChemSpID
	childrenProc<-sapply(1:length(CTS.options), function(i)
	{
		cat(CTS.options[i],"\n")
		CTSgetR(id=id,to=CTS.options[i],from="ChemSpider")
	})
	print(childrenProc)
	
	require(RJSON)
	urlInchi <- paste("http://cts.fiehnlab.ucdavis.edu/service/compound/", key, sep='')
	JSONstring <- getURL(urlInchi)
	Content <- fromJSON(JSONstring, method = "C", unexpected.escape = "error" )
	#print(Content)
  # Postprocess:
  # Split CAS, SID, CID, IUPAC
  # (don't split names yet, since we don't have a good rule. - and , 
  # are both problematic here)
 
  tosplit <- c("cas", "sid", "cid", "iupac", "smiles", "kegg", "chebi")
  #for(var in tosplit)
  #{
  #  if(childrenProc[[var]] != "error")
  #  {
  #    childrenProc[[var]] <- strsplit(childrenProc[[var]], ", ", fixed=TRUE)
  #    childrenProc[[var]] <- unlist(lapply(childrenProc[[var]],
  #      function(x) sub("^ *([^ ]+) *$", "\\1", x)))
  #  } 
  #}
  # Try to handle names in a graceful way
  #if(!is.null(childrenProc$names))
  #{
    # add final comma and space to match the last name correctly
  #  names <- paste(childrenProc$names, ", ", sep='')
  #  matches.list <- gregexpr(" (.*?) - (-*[0-9]+), ", names, perl=TRUE)
  #  matches <- regmatches(names, matches.list )[[1]]
  #  matches <- t(sapply(matches, function(match)
  #  {
  #    match.name <- sub(" (.*?) - (-*[0-9]+), ", "\\1", match)
  #    match.score <- as.integer(sub(" (.*?) - (-*[0-9]+), ", "\\2", match))
  #   return(list(name=match.name, score=match.score))
  #  }, USE.NAMES=FALSE))
  #  childrenProc$names <- matches
    # " (.*?) - (-*[0-9]+), "
  #}

  # Check and fix CAS (eliminate the 12-34-5CHEBI and NIKKAJI stuff)
  if(!is.null(childrenProc$CAS))
  {
    childrenProc$CAS <- childrenProc$CAS[which(grepl("^[-0-9]+$", childrenProc$CAS))]
  }
  print(Content)
  return(Content)
}