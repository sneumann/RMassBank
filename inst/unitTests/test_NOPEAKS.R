test.nopeaks <- function(){
	w <- newMsmsWorkspace()
	w@aggregated <- data.frame(mzFound = numeric(0), intensity = numeric(0), good = logical(0), mzCalc = numeric(0), formula = character(0), dbe = numeric(0),
		formulaCount = integer(0), dppm = numeric(0), dppmBest = numeric(0), scan = integer(0), cpdID = character(0), parentScan = integer(0), dppmRc = numeric(0),
		index = integer(0), noise = logical(0), reanalyzed.formula = character(0), reanalyzed.mzCalc = numeric(0), reanalyzed.dppm = numeric(0),reanalyzed.formulaCount = numeric(0),
		reanalyzed.dbe = numeric(0),matchedReanalysis = logical(0))
	filterMultiplicity(w, mode="pH")
	checkTrue(TRUE)
}