library(RMassBank)
file.copy(system.file("list/NarcoticsDataset.csv",package="RMassBankData"), "./Compoundlist.csv")

loadRmbSettings(system.file("XCMSinput/mysettings.ini",package="RMassBank"))
w <- newMsmsWorkspace()
files <- list.files(system.file("spectra", package="RMassBankData"), ".mzML", full.names = TRUE)
w@files <- files[1]
loadList("./Compoundlist.csv")
w <- msmsWorkflow(w, mode="pH", steps=c(1:8))