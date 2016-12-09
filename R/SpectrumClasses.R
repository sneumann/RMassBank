#' @import MSnbase
#' @import Biobase
#' @import S4Vectors

#' @exportClass RmbSpectrum2
.RmbSpectrum2 <- setClass("RmbSpectrum2",
		representation = representation(
				satellite="logical",
				low="logical",
				rawOK ="logical",
				good = "logical",
				mzCalc = "numeric",
				formula = "character",
				dbe = "numeric",
				formulaCount = "integer",
				formulaSource = "character",
				dppm = "numeric",
				dppmBest = "numeric",
				ok = "logical",
				info = "list",
				properties = "data.frame"
		),
		contains=c("Spectrum2"),
		prototype = prototype(
				satellite = logical(),
				low = logical(),
				rawOK = logical(),
				good = logical(),
				mzCalc = numeric(),
				formula = character(),
				dbe = numeric(),
				formulaCount = integer(),
				formulaSource = character(),
				dppm = numeric(),
				dppmBest = numeric(),
				ok = logical(),
				info = list(),
				properties = data.frame(),
				new("Versioned", versions=c(classVersion("Spectrum2"), RmbSpectrum2 = "0.1.2"))
		),
)

#' @exportClass RmbSpectrum2List
.RmbSpectrum2List <- setClass("RmbSpectrum2List", contains="SimpleList",
		prototype=prototype(elementType="RmbSpectrum2"))
#
#setAs("ANY", "RmbSpectrum2List", function(from) {
#			coerceToSimpleList(from)
#		})

#' @exportClass RmbSpectraSet
.RmbSpectraSet <- setClass("RmbSpectraSet",
		representation = representation(
				parent = "Spectrum1",
				children = "RmbSpectrum2List",
				# These are done as slots and not as S4 functions, because they are set during the workflow
				# in "checking" steps. It's easier.
				found = "logical",
				complete = "logical",
				empty = "logical",
				formula = "character",
				id = "character",
				mz = "numeric",
				name = "character",
				mode = "character",
        smiles = "character"
				#annotations = "list"
				),
    contains=c("Versioned"),
		prototype = prototype(
				parent = new("Spectrum1"),
				children = new("RmbSpectrum2List"),
				found = FALSE,
				complete = NA,
				empty = NA,
				formula = character(),
				id = character(),
				mz = numeric(),
				name = character(),
				mode = character(),
        smiles = character(),
        new("Versioned", versions=c(RmbSpectraSet = "0.1.2"))
    # version 0.1.1: introduced versioning and SMILES slot
    # version 0.1.2: feed polarity to parent and children
		)
);

#' @exportClass RmbSpectraSetList
.RmbSpectraSetList <- setClass("RmbSpectraSetList", contains="SimpleList",
		prototype=prototype(elementType="RmbSpectraSet"))




setGeneric("getData",	function(s) standardGeneric("getData"))
setGeneric("setData",	function(s, df, ...) standardGeneric("setData"))

