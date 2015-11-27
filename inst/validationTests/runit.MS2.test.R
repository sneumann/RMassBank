test.NA <- function(){
	RUnit::checkTrue(!(NA %in% as.matrix(RMassBank.env$mb[[RMassBank.env$testnumber]]@compiled_ok[[1]][['PK$PEAK']])))
}

test.peaksvsprecursor <- function(){
	Max_Peak <- unname(RMassBank.env$mb[[RMassBank.env$testnumber]]@compiled_ok[[1]][['PK$PEAK']][dim(RMassBank.env$mb[[RMassBank.env$testnumber]]@compiled_ok[[1]][['PK$PEAK']])[1],1])
	Precursor <- RMassBank.env$mb[[RMassBank.env$testnumber]]@compiled_ok[[1]][['MS$FOCUSED_ION']][['PRECURSOR_M/Z']]
	if(is.na(Precursor)){
		RUnit::checkTrue(TRUE)
	}else{
		RUnit::checkEquals(Max_Peak, Precursor, tolerance = Precursor/100)
	}
}

test.precursormz <- function(){
	precursorlist <- c("[M+H]+","[M+Na]+","[M-H]-","[M+HCOO-]-","[M]+","[M]-")
	if(is.na(RMassBank.env$mb[[RMassBank.env$testnumber]]@compiled_ok[[1]][['MS$FOCUSED_ION']][['PRECURSOR_TYPE']]) || is.na(RMassBank.env$mb[[RMassBank.env$testnumber]]@compiled_ok[[1]][['MS$FOCUSED_ION']][['PRECURSOR_M/Z']])){
		RUnit::checkTrue(TRUE)
	} else{
		precursor <- grep(RMassBank.env$mb[[RMassBank.env$testnumber]]@compiled_ok[[1]][['MS$FOCUSED_ION']][['PRECURSOR_TYPE']],precursorlist, value = TRUE, fixed = TRUE)
		if(precursor == "[M+H]+"){
		RUnit::checkEquals(RMassBank.env$mb[[RMassBank.env$testnumber]]@compiled_ok[[1]][['MS$FOCUSED_ION']][['PRECURSOR_M/Z']],RMassBank.env$mb[[RMassBank.env$testnumber]]@compiled_ok[[1]][['CH$EXACT_MASS']] + 1.008,tolerance = 0.002)
		}
		if(precursor == "[M+Na]+"){
			RUnit::checkEquals(RMassBank.env$mb[[RMassBank.env$testnumber]]@compiled_ok[[1]][['MS$FOCUSED_ION']][['PRECURSOR_M/Z']],RMassBank.env$mb[[RMassBank.env$testnumber]]@compiled_ok[[1]][['CH$EXACT_MASS']] + 22.989,tolerance = 0.002)
		}
		if(precursor == "[M-H]-"){
			RUnit::checkEquals(RMassBank.env$mb[[RMassBank.env$testnumber]]@compiled_ok[[1]][['MS$FOCUSED_ION']][['PRECURSOR_M/Z']],RMassBank.env$mb[[RMassBank.env$testnumber]]@compiled_ok[[1]][['CH$EXACT_MASS']] - 1.008,tolerance = 0.002)
		}
		if(precursor == "[M+HCOO-]-"){
			RUnit::checkEquals(RMassBank.env$mb[[RMassBank.env$testnumber]]@compiled_ok[[1]][['MS$FOCUSED_ION']][['PRECURSOR_M/Z']],RMassBank.env$mb[[RMassBank.env$testnumber]]@compiled_ok[[1]][['CH$EXACT_MASS']] + 45.017,tolerance = 0.002)
		}
	}
}

test.PrecursorType <- function(){    
	precursorlist <- c("[M+H]+","[M+Na]+","[M-H]-","[M+HCOO-]-","[M]+","[M]-")
	if(is.na(RMassBank.env$mb[[RMassBank.env$testnumber]]@compiled_ok[[1]][['MS$FOCUSED_ION']][['PRECURSOR_TYPE']])){
		RUnit::checkTrue(TRUE)
	}else{
	RUnit::checkTrue(RMassBank.env$mb[[RMassBank.env$testnumber]]@compiled_ok[[1]][['MS$FOCUSED_ION']][['PRECURSOR_TYPE']] %in% precursorlist)
	}
}

test.smilesvsexactmass <- function(){
	Mass_Calculated_Through_Smiles <- smiles2mass(RMassBank.env$mb[[RMassBank.env$testnumber]]@compiled_ok[[1]][['CH$SMILES']])
	Exact_Mass <- RMassBank.env$mb[[RMassBank.env$testnumber]]@compiled_ok[[1]][['CH$EXACT_MASS']]
	RUnit::checkEquals(Mass_Calculated_Through_Smiles, Exact_Mass, tolerance = Exact_Mass/100)
}

test.sumintensities <- function(){
	sumOfIntensities <- sum(RMassBank.env$mb[[RMassBank.env$testnumber]]@compiled_ok[[1]][['PK$PEAK']][,2])
	RUnit::checkTrue(sumOfIntensities > 0)
}

test.TitleVsType <- function(){
		RMassBank.env$testnumber <- RMassBank.env$testnumber + 1
		if(is.na(RMassBank.env$mb[[RMassBank.env$testnumber-1]]@compiled_ok[[1]][['MS$FOCUSED_ION']][['PRECURSOR_TYPE']])){
			RUnit::checkTrue(TRUE)
		}else{
			RUnit::checkTrue(grepl(RMassBank.env$mb[[RMassBank.env$testnumber-1]]@compiled_ok[[1]][['MS$FOCUSED_ION']][['PRECURSOR_TYPE']], RMassBank.env$mb[[RMassBank.env$testnumber-1]]@compiled_ok[[1]][['RECORD_TITLE']], fixed = TRUE))
		}
}
