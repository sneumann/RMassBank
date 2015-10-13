test.mzRRead <- function(){
	allOK <- TRUE
	records <- list.files("XX/recdata",full.names=TRUE)
	rightrecords <- list.files(system.file("records/XX/recdata", package="RMassBankData"), full.names=TRUE)
	
	for(i in 1:length(records)){
		gen <- file(description = records[i], open = "r", blocking = TRUE,
		encoding = getOption("encoding"), raw = FALSE)
		genLines <- readLines(gen)
		RMBLine <- grep("MS$DATA_PROCESSING: WHOLE",genLines, fixed=TRUE)
        DATELine <- grep("DATE: ",genLines, fixed=TRUE)
		genLines <- genLines[-c(RMBLine,DATELine)]
		print(rightrecords[i])
		rig <- file(description = rightrecords[i], open = "r", blocking = TRUE,
		encoding = getOption("encoding"), raw = FALSE)
		rigLines <- readLines(rig)
		RMBLine <- grep("MS$DATA_PROCESSING: WHOLE",rigLines, fixed=TRUE)
        DATELine <- grep("DATE: ",rigLines, fixed=TRUE)
		rigLines <- rigLines[-c(RMBLine,DATELine)]
		close(gen)
		close(rig)
		if(!identical(genLines,rigLines)){
			print(i)
			allOK <- FALSE
			break
		}

	}
	RUnit::checkTrue(allOK)
}