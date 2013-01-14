test.instrumentname <- function(){
	Instrument_Name <- mb@compiled_ok[[testNumber]][['AC$INSTRUMENT']]
	checkTrue(Instrument_Name %in% Instrument_List)
}

test.NA <- function(){
	checkTrue(!(NA %in% as.matrix(mb@compiled_ok[[testNumber]][['PK$PEAK']])))
}

test.peaksvsprecursor <- function(){
	Max_Peak <- unname(mb@compiled_ok[[testNumber]][['PK$PEAK']][dim(mb@compiled_ok[[testNumber]][['PK$PEAK']])[1],1])
	Precursor <- mb@compiled_ok[[testNumber]][['MS$FOCUSED_ION']][['PRECURSOR_M/Z']]
	if(is.na(Precursor)){
		checkTrue(TRUE)
	}else{
		checkEquals(Max_Peak, Precursor, tolerance = Precursor/100)
	}
}

test.precursormz <- function(){
	precursorlist <- c("[M+H]+","[M+Na]+","[M-H]-","[M+HCOO-]-","[M]+","[M]-")
	if(is.na(mb@compiled_ok[[testNumber]][['MS$FOCUSED_ION']][['PRECURSOR_TYPE']])){
		checkTrue(TRUE)
	} else{
		precursor <- grep(mb@compiled_ok[[testNumber]][['MS$FOCUSED_ION']][['PRECURSOR_TYPE']],precursorlist, value = TRUE, fixed = TRUE)
		if(precursor == "[M+H]+"){
		checkEquals(mb@compiled_ok[[testNumber]][['MS$FOCUSED_ION']][['PRECURSOR_M/Z']],mb@compiled_ok[[testNumber]][['CH$EXACT_MASS']] + 1.008,tolerance = 0.002)
		}
		if(precursor == "[M+Na]+"){
			checkEquals(mb@compiled_ok[[testNumber]][['MS$FOCUSED_ION']][['PRECURSOR_M/Z']],mb@compiled_ok[[testNumber]][['CH$EXACT_MASS']] + 22.989,tolerance = 0.002)
		}
		if(precursor == "[M-H]-"){
			checkEquals(mb@compiled_ok[[testNumber]][['MS$FOCUSED_ION']][['PRECURSOR_M/Z']],mb@compiled_ok[[testNumber]][['CH$EXACT_MASS']] - 1.008,tolerance = 0.002)
		}
		if(precursor == "[M+HCOO-]-"){
			checkEquals(mb@compiled_ok[[testNumber]][['MS$FOCUSED_ION']][['PRECURSOR_M/Z']],mb@compiled_ok[[testNumber]][['CH$EXACT_MASS']] + 45.017,tolerance = 0.002)
		}
	}
}

test.PrecursorType <- function(){    
	precursorlist <- c("[M+H]+","[M+Na]+","[M-H]-","[M+HCOO-]-","[M]+","[M]-")
	if(is.na(mb@compiled_ok[[testNumber]][['MS$FOCUSED_ION']][['PRECURSOR_TYPE']])){
		checkTrue(TRUE)
	}else{
	checkTrue(mb@compiled_ok[[testNumber]][['MS$FOCUSED_ION']][['PRECURSOR_TYPE']] %in% precursorlist)
	}
}

test.smilesvsexactmass <- function(){
	Mass_Calculated_Through_Smiles <- smiles2mass(mb@compiled_ok[[testNumber]][['CH$SMILES']])
	Exact_Mass <- mb@compiled_ok[[testNumber]][['CH$EXACT_MASS']]
	checkEquals(Mass_Calculated_Through_Smiles, Exact_Mass, Exact_Mass/100)
}

test.sumintensities <- function(){
	sumOfIntensities <- sum(mb@compiled_ok[[testNumber]][['PK$PEAK']][,2])
	checkTrue(sumOfIntensities > 0)
}

test.TitleVsType <- function(){
		testNumber <<- testNumber + 1
		if(is.na(mb@compiled_ok[[testNumber-1]][['MS$FOCUSED_ION']][['PRECURSOR_TYPE']])){
			checkTrue(TRUE)
		}else{
		checkTrue(grepl(mb@compiled_ok[[testNumber-1]][['MS$FOCUSED_ION']][['PRECURSOR_TYPE']], mb@compiled_ok[[testNumber-1]][['RECORD_TITLE']], fixed = TRUE))
		}
}
