setClass("RmbSpectraSet",
		representation = representation(
				parent = "Spectrum1",
				children = "list",
				# These are done as slots and not as S4 functions, because they are set during the workflow
				# in "checking" steps. It's easier.
				found = "logical",
				complete = "logical",
				empty = "logical",
				formula = "character",
				id = "integer",
				mz = "numeric",
				name = "character",
				annotations = "list"
				),
		prototype = prototype(
				parent = new("Spectrum1"),
				children = list(),
				found = FALSE,
				complete = NA,
				empty = NA,
				formula = character(),
				id = integer(),
				mz = numeric(),
				name = character(),
				annotations = list()
		)
		,
		validity = function(object)
		{
			childrenSpectraOK <- all(unlist(lapply(object@children, function(s) inherits(s, "RmbSpectrum2"))))
			if(!childrenSpectraOK) return("MS2 spectra are not of class RmbSpectrum2")
			return(TRUE)
		}
);

setClass("RmbSpectrum2",
		representation = representation(
				satellite="logical",
				low="logical",
				rawOK ="logical",
				good = "logical",
				mzCalc = "numeric",
				formula = "numeric",
				dppm = "numeric"
				),
		contains=c("Spectrum2"),
		prototype = prototype(
				satellite = logical(),
				low = logical(),
				rawOK = logical(),
				good = logical(),
				mzCalc = numeric(),
				formula = numeric(),
				dppm = numeric()
				),
)

setGeneric("getData",	function(s) standardGeneric("getData"))
setGeneric("setData",	function(s, df) standardGeneric("setData"))

