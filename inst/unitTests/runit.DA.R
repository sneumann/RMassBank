test.mzRRead <- function(){
	allOK <- TRUE
	records <- list.files("XX/recdata",full.names=TRUE)
	rightrecords <- list.files(system.file("records/XX/recdata/", package="RMassBankData"), full.names=TRUE)
	
	for(i in 1:length(records)){
		gen <- file(description = records[i], open = "r", blocking = TRUE,
		encoding = getOption("encoding"), raw = FALSE)
		genLines <- readLines(gen)
		RMBLine <- grep("MS$DATA_PROCESSING: WHOLE",genLines, fixed=TRUE)
		genLines <- genLines[-RMBLine]
		rig <- file(description = rightrecords[i], open = "r", blocking = TRUE,
		encoding = getOption("encoding"), raw = FALSE)
		rigLines <- readLines(rig)
		RMBLine <- grep("MS$DATA_PROCESSING: WHOLE",rigLines, fixed=TRUE)
		rigLines <- rigLines[-RMBLine]
		close(gen)
		close(rig)
		if(!identical(genLines,rigLines)){
			allOK <- FALSE
			break
		}

	}
	checkTrue(allOK)
}