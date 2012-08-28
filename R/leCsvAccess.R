
.listEnvEnv <- new.env()
assign("listEnv", NULL, envir=.listEnvEnv)

#' Load compound list for RMassBank
#' 
#' Loads a CSV compound list with compound IDs
#' 
#' The list is loaded into the variable \code{\var{compoundList}} in the environment
#' \code{listEnv} (which defaults to the global environment) and used by
#' the \code{findMz}, \code{findCAS}, ... functions.
#' 
#' resetList() clears a currently loaded list.
#' 
#' @aliases loadList resetList
#' @usage loadList(path, listEnv=NULL)
#' 
#' 			resetList()
#' @param path Path to the CSV list.
#' @param listEnv The environment to load the list into. By default, the namelist is loaded
#' 		into an environment internally in RMassBank. 
#' @return No return value.
#' @author Michael Stravs
#' @seealso \code{\link{findMz}}
#' @examples
#' 
#' ##
#' \dontrun{loadList("mylist.csv")}
#' 
#' @export
loadList <- function(path, listEnv = NULL)
{
	if(is.null(listEnv))
		listEnv <- .listEnvEnv
  compoundList <- read.csv(path, stringsAsFactors=FALSE)
  # check whether the necessary columns are present
  n <- colnames(compoundList)
  cols <- c('ID', 'Name', 'SMILES', 'RT', 'CAS')
  d <- setdiff(cols, n)
  if(length(d)>0)
      stop("Some columns are missing in the compound list. Needs at least ID, Name, SMILES, RT.")
  if(length(duplicated(compoundList$cpdID)) > 0)
      stop("Duplicate compound IDs are present. Please check your list.")
  assign("listEnv", listEnv, envir=.listEnvEnv) 
  .listEnvEnv$listEnv$compoundList <- compoundList
}

#' @export
resetList <- function()
{
	if(is.null(.listEnvEnv$listEnv))
		return()
	if(!exists("compoundList", where=.listEnvEnv$listEnv))
		return()
	rm(.listEnvEnv$listEnv$compoundList)
	assign("listEnv", NULL, envir=.listEnvEnv)
}

#' Create Rcdk molecule from SMILES
#' 
#' Generates a Rcdk molecule object from SMILES code, which is fully typed and
#' usable (in contrast to the built-in \code{parse.smiles}).
#' 
#' \bold{NOTE: As of today (2012-03-16), Rcdk discards stereochemistry when
#' loading the SMILES code!} Therefore, do not trust this function blindly,
#' e.g.  don't generate InChI keys from the result. It is, however, useful if
#' you want to compute the mass (or something else) with Rcdk.
#' 
#' @usage getMolecule(smiles)
#' @param smiles The SMILES code of the compound.
#' @return A Rcdk \code{IAtomContainer} reference.
#' @author Michael Stravs
#' @seealso \code{\link{parse.smiles}}
#' @examples
#' 
#' # Lindane:
#' getMolecule("C1(C(C(C(C(C1Cl)Cl)Cl)Cl)Cl)Cl")
#' # Benzene:
#' getMolecule("C1=CC=CC=C1")
#' 
#' @export
getMolecule <- function(smiles)
{
  capture.output(mol <- parse.smiles(smiles)[[1]])
  do.aromaticity(mol)
  convert.implicit.to.explicit(mol)
  do.aromaticity(mol)
  do.typing(mol)
  do.isotopes(mol)
 
  return(mol)
}



#' Find the exact mass +/- a given margin for a given formula or its ions and adducts.
#' 
#' @param formula The molecular formula  in text or list format (see \code{\link{formulastring.to.list}}
#' @param mode \code{"pH", "pNa", "pM", "mH", "mM", "mFA"} for different ions 
#' 			([M+H]+, [M+Na]+, [M]+, [M-H]-, [M]-, [M+FA]-). "" for the uncharged molecule.
#' @param ppm The ppm margin to add/subtract
#' @param deltaMz The absolute mass to add/subtract. Cumulative with \code{ppm}
#' @return A \code{list(mzMin=, mzCenter=, mzMax=)} with the masses.
#' @author Michael A. Stravs, Eawag <michael.stravs@@eawag.ch>
#' @examples findMz.formula("C6H6")
#' @seealso \code{\link{findMz}}
#' @export
findMz.formula <- function(formula, mode="pH", ppm=10, deltaMz=0)
{
	mzopt <- list(addition="", charge=0)
	if(mode=="pH") mzopt <- list(addition="H", charge=1)
	if(mode=="pNa") mzopt <- list(addition="Na", charge=1)
	if(mode=="pM") mzopt <- list(addition="", charge=1)
	if(mode=="mH") mzopt <- list(addition="H-1", charge=-1)
	if(mode=="mFA") mzopt <- list(addition="C1O2", charge=-1)
	if(mode=="mM") mzopt <- list(addition="", charge=-1)
	if(mode=="") mzopt <- list(addition="", charge=0)
	
	formula <- add.formula(formula, mzopt$addition)
	formula <- get.formula(formula, charge=mzopt$charge)
	m <- formula@mass
	delta <- ppm(m, ppm, l=TRUE)
	return(list(mzMin=delta[[2]] - deltaMz, mzMax=delta[[1]] + deltaMz, mzCenter=m))
}

#' Find compound information
#' 
#' Retrieves compound information from the loaded compound list or calculates
#' it from the SMILES code in the list.
#' 
#' @aliases findMz findSmiles findFormula findRt findCAS findName
#' @usage  findMz(cpdID, mode = "pH", ppm = 10, deltaMz = 0)
#' 
#' findRt(cpdID) 
#' 
#' findSmiles(cpdID) 
#' 
#' findFormula(cpdID) 
#' 
#' findCAS(cpdID)
#' 
#' findName(cpdID)
#' @param cpdID The compound ID in the compound list.
#' @param mode Specifies the species of the molecule: An empty string specifies
#' uncharged monoisotopic mass, \code{\var{pH}} (positive H) specifies [M+H]+,
#' \code{\var{pNa}} specifies [M+Na]+, \code{\var{pM}} specifies [M]+,
#' \code{\var{mH}} and \code{\var{mFA}} specify [M-H]- and [M+FA]-,
#' respectively. (I apologize for the naming of \code{\var{pH}} which has
#' absolutely nothing to do with chemical \emph{pH} values.)
#' @param ppm Specifies ppm window (10 ppm will return the range of the
#' molecular mass + and - 10 ppm).
#' @param deltaMz Specifies additional m/z window to add to the range (deltaMz
#' = 0.02 will return the range of the molecular mass +- 0.02 (and additionally
#' +- the set ppm value).
#' @return \code{findMz} will return a \code{list(mzCenter=, mzMin=, mzMax=)}
#' with the molecular weight of the given ion, as calculated from the SMILES
#' code and Rcdk.
#' 
#' \code{findRt}, \code{findSmiles},\code{findCAS},\code{findName} will return
#' the corresponding entry from the compound list. \code{findFormula} returns
#' the molecular formula as determined from the SMILES code.
#' 
#' @author Michael Stravs
#' @seealso \code{\link{findMass}}, \code{\link{loadList}}, \code{\link{findMz.formula}}
#' @examples
#' 
#' \dontrun{%
#' 	findMz(123, "pH", 5)
#' 	findFormula(123)
#' }
#' 
#' @export
findMz <- function(cpdID, mode="pH", ppm=10, deltaMz=0)
{
  s <- findSmiles(cpdID)
  if(s=="-") s <- NA
  if(is.na(s)) return(list(mzMin=NA,mzMax=NA,mzCenter=NA))
  formula <- .get.mol2formula(getMolecule(s))
  return(findMz.formula(formula@string, mode, ppm, deltaMz))
}

#findMz <- function(...)	findMz.rcdk(...)
#' @export
findRt <- function(cpdID) {
	if(is.null(.listEnvEnv$listEnv))
		stop("Compound list must be loaded first.")
	if(!exists("compoundList", where=.listEnvEnv$listEnv))
		stop("Compound list must be loaded first.")
	if(is.character(cpdID))
		cpdID <- as.numeric(cpdID)
	rt <- .listEnvEnv$listEnv$compoundList[which(.listEnvEnv$listEnv$compoundList$ID == cpdID),"RT"]
	if(!is.null(getOption("RMassBank")$rtShift))
		rt <- rt + getOption("RMassBank")$rtShift
	if(is.na(rt)) return(list(RT=NA))
	else return(list(RT=rt))
}

#' @export
findSmiles <- function(cpdID) {
	if(is.character(cpdID))
		cpdID <- as.numeric(cpdID)
	if(is.null(.listEnvEnv$listEnv))
		stop("Compound list must be loaded first.")
	if(!exists("compoundList", where=.listEnvEnv$listEnv))
		stop("Compound list must be loaded first.")
	
	return(.listEnvEnv$listEnv$compoundList[match(cpdID, .listEnvEnv$listEnv$compoundList$ID),"SMILES"])
}

#' @export
findFormula <- function(cpdID) {
  smiles <- findSmiles(cpdID)
  mol <- getMolecule(smiles)
  f <- .get.mol2formula(mol)
  return(f@string)
}

#' @export
findCAS <- function(cpdID) {
	if(is.character(cpdID))
		cpdID <- as.numeric(cpdID)
	if(is.null(.listEnvEnv$listEnv))
		stop("Compound list must be loaded first.")
	if(!exists("compoundList", where=.listEnvEnv$listEnv))
		stop("Compound list must be loaded first.")
	return(.listEnvEnv$listEnv$compoundList[which(.listEnvEnv$listEnv$compoundList$ID == cpdID),"CAS"])
}

#' @export
findName <- function(cpdID) {
	if(is.character(cpdID))
		cpdID <- as.numeric(cpdID)
	if(is.null(.listEnvEnv$listEnv))
		stop("Compound list must be loaded first.")
	if(!exists("compoundList", where=.listEnvEnv$listEnv))
		stop("Compound list must be loaded first.")
	return(.listEnvEnv$listEnv$compoundList[which(.listEnvEnv$listEnv$compoundList$ID == cpdID),"Name"])
}

#' Calculate exact mass
#' 
#' Retrieves the exact mass of the uncharged molecule. It works directly from
#' the SMILES and therefore is used in the MassBank workflow
#' (\code{\link{mbWorkflow}}) - there, all properties are calculated from the
#' SMILES code retrieved from the database. (Alternatively, takes also the
#' compound ID as parameter and looks it up.) Calculation relies on Rcdk.
#' 
#' @usage findMass(cpdID_or_smiles)
#' @param cpdID_or_smiles SMILES code or compound ID of the molecule. (Numerics
#' are treated as compound ID).
#' @return  Returns the exact mass of the uncharged molecule.
#' @author Michael Stravs
#' @seealso \code{\link{findMz}}
#' @examples
#' 
#' ##
#' findMass("OC[C@@H]1OC(O)[C@@H](O)[C@@@@H](O)[C@@@@H]1O")
#' 
#' @export
findMass <- function(cpdID_or_smiles)
{
	if(!is.numeric(cpdID_or_smiles))
		s <- cpdID_or_smiles
	else
		s <- findSmiles(cpdID_or_smiles)
	mol <- getMolecule(s)
	return(get.exact.mass(mol))
}

