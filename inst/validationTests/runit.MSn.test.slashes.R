test.slashes <- function(){
	Type <- as.numeric(substring(mb@compiled_ok[[testNumber]][['AC$MASS_SPECTROMETRY']][['MS_TYPE']],first = 3))
	slashes <- length(gregexpr('/', mb@compiled_ok[[testNumber]][['MS$FOCUSED_ION']][['PRECURSOR_M/Z']]))
	RUnit::checkEquals(Type - 2,slashes)
}