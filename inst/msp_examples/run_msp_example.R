library(RMassBank)
w <- newMsmsWorkspace()
files <- list.files(system.file('msp_examples', package="RMassBankData"), '.msp', full.names=TRUE)
w@files <- files
loadList(system.file('msp_examples/Compoundlist.csv', package="RMassBankData"))
loadRmbSettings(system.file('msp_examples/RMB_options.ini', package="RMassBankData"))
w <- msmsWorkflow(w, readMethod='msp', 
                  filetable=system.file('msp_examples/Filelist.csv', package="RMassBankData"),
                  mode='pH', steps=1, archivename='msp_archive')
mb <- newMbWorkspace(w)
#mb <- mbWorkflow(mb)
mb <- resetInfolists(mb)
mb <- loadInfolists(mb, system.file('infolists', package="RMassBankData"))
mb <- mbWorkflow(mb, filter=FALSE)
