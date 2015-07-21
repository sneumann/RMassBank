#' @import MSnbase
#' @import Biobase
#' @import S4Vectors
#' @importClassesFrom MSnbase Spectrum Spectrum2 


# @exportClass RmbSpectrum2
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
				dppm = "numeric",
				dppmBest = "numeric",
				ok = "logical",
				info = "list"
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
				dppm = numeric(),
				dppmBest = numeric(),
				ok = logical(),
				info = list(),
				new("Versioned", versions=c(classVersion("Spectrum2"), RmbSpectrum2 = "0.1.0"))
		),
)

# @exportClass RmbSpectrum2List
.RmbSpectrum2List <- setClass("RmbSpectrum2List", contains="SimpleList",
		prototype=prototype(elementType="RmbSpectrum2"))
#
#setAs("ANY", "RmbSpectrum2List", function(from) {
#			coerceToSimpleList(from)
#		})

# @exportClass RmbSpectraSet
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
				mode = "character"
				#annotations = "list"
				),
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
				mode = character()
		)
);

# @exportClass RmbSpectraSetList
.RmbSpectraSetList <- setClass("RmbSpectraSetList", contains="SimpleList",
		prototype=prototype(elementType="RmbSpectraSet"))




setGeneric("getData",	function(s) standardGeneric("getData"))
setGeneric("setData",	function(s, df, ...) standardGeneric("setData"))

