test.eletronicnoise <- function(){
	failpeaks <- read.csv("pH_narcotics_Failpeaks.csv")
	RUnit::checkTrue(!any(apply(failpeaks,1,function(x) all(x[2:5] == c(1738,2819,668,201.69144)))))
}

test.formulacalculation <- function(){
	failpeaks <- read.csv("pH_narcotics_Failpeaks.csv")
	
	RUnit::checkTrue(!any(apply(failpeaks,1,function(x) all(x[2:5] == c(70,2758,321,56.04933)))))
}