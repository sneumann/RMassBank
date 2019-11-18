# TODO: Add comment
# 
# Author: stravsmi
###############################################################################


# TODO: Add comment

#' Extract an MS/MS spectrum from MS2 TIC
#' 
#' Extract an MS/MS spectrum or multiple MS/MS spectra based on the TIC of the MS2 and precursor mass, picking
#' the most intense MS2 scan. Can be used, for example, to get a suitable MS2 from direct infusion data which was
#' collected with purely targeted MS2 without MS1.
#' 
#' Note that this is not a precise function and only really makes sense in direct infusion and if
#' the precursor is really known, because MS2 precursor data is only "roughly" accurate (to 2 dp).
#' The regular \code{findMsMsHR} functions confirm the exact mass of the precursor in the MS1 scan.
#' 
#' @aliases findMsMsHR.ticMS2
#' @param msRaw The mzR raw file
#' @param mz Mass to find
#' @param limit.coarse Allowed mass deviation for scan precursor (in m/z values)
#' @param limit.fine Unused here, but present for interface compatiblity with findMsMsHR
#' @param rtLimits Unused here, but present for interface compatiblity with findMsMsHR
#' @param maxCount Maximal number of spectra to return
#' @param headerCache Cached results of header(msRaw), either to speed up the operations or to operate with
#' 			preselected header() data 
#' @param fillPrecursorScan Unused here, but present for interface compatiblity with findMsMsHR
#' @param deprofile Whether deprofiling should take place, and what method should be
#' 			used (cf. \code{\link{deprofile}}) 
#' @param trace Either \code{"ms2tic"} or \code{"ms2basepeak"}: Which intensity trace to use - can be either the
#' 			TIC of the MS2 or the basepeak intensity of the MS2.
#' @return	a list of "spectrum sets" as defined in \code{\link{findMsMsHR}}, sorted
#' 			by decreasing precursor intensity.
#' 
#' @author stravsmi
#' @export
findMsMsHR.ticms2 <- function(msRaw, mz, limit.coarse, limit.fine, rtLimits = NA, maxCount = NA,
    headerCache = NULL, fillPrecursorScan = FALSE,
    deprofile = getOption("RMassBank")$deprofile, trace = "ms2tic")
{
  #	if(!is.na(rtLimits))
  #	{  
  #		eic <- subset(eic, rt >= rtLimits[[1]] & rt <= rtLimits[[2]])
  #	}
  if(!is.null(headerCache))
    headerData <- headerCache
  else
    headerData <- as.data.frame(header(msRaw))
  
  # find MS2 tic or MS2 base peak chromatogram
  #return(data.frame(rt = rt, intensity = pks_t, scan = scan))
  if(trace == "ms2tic")
    eic <- headerData[,c("retentionTime", "totIonCurrent", "acquisitionNum")] 
  else if(trace == "ms2basepeak")
    eic <- headerData[,c("retentionTime", "basePeakIntensity", "acquisitionNum")]
  else
    stop("Select a valid trace parameter (ms2tic or ms2basepeak).")
  colnames(eic) <- c("rt", "intensity", "scan")
  
  
  # Find MS2 spectra with precursors which are in the allowed 
  # scan filter (coarse limit) range
  findValidPrecursors <- headerData[
      which(headerData$precursorMZ > (mz - limit.coarse) &
            headerData$precursorMZ < (mz + limit.coarse)),]
  # Find the precursors for the found spectra
  
  
  # Crop the "EIC" to the valid scans
  eic <- eic[eic$scan %in% findValidPrecursors$acquisitionNum,]
  # Order by intensity, descending
  eic <- eic[order(eic$intensity, decreasing=TRUE),]
  if(nrow(eic) == 0)
    return(list(list(foundOK = FALSE)))
  if(!is.na(maxCount))
  {
    spectraCount <- min(maxCount, nrow(eic))
    eic <- eic[1:spectraCount,]
  }
  # Construct all spectra groups in decreasing intensity order
  spectra <- lapply(eic$scan, function(masterScan)
      {
        ret <- list()
        ret$parentHeader <- matrix(0, ncol = 20, nrow = 1)
        
        rownames(ret$parentHeader) <- 1
        colnames(ret$parentHeader) <- c("seqNum", "acquisitionNum", "msLevel", "peaksCount", "totIonCurrent", "retentionTime", "basepeakMZ", 
            "basePeakIntensity", "collisionEnergy", "ionisationEnergy", "lowMZ", "highMZ", "precursorScanNum",
            "precursorMZ", "precursorCharge", "precursorIntensity", "mergedScan", "mergedResultScanNum", 
            "mergedResultStartScanNum", "mergedResultEndScanNum")
        ret$parentHeader[1,1:3] <- 1
        ##Write nothing in the parents if there is no MS1-spec
        ret$parentHeader[1,4:20] <- 0
        ret$parentHeader[1,6] <- NA

        childHeaders <- headerData[
            which(headerData$acquisitionNum == masterScan 
                  & headerData$precursorMZ > (mz - limit.coarse) 
                  & headerData$precursorMZ < (mz + limit.coarse)), ,
            drop = FALSE]
      
        childScans <- childHeaders$acquisitionNum
        ret$parentScan <- min(childScans)-1
		ret$parentHeader[1,1:3] <- min(childScans)-1
		
        ret$parentPeak <- matrix(nrow = 1, ncol = 2)
        colnames(ret$parentPeak) <- c("mz","int")
        ret$parentPeak[1,] <- c(mz,100)
        
        msmsPeaks <- lapply(childHeaders$seqNum, function(scan)
            {
              pks <- mzR::peaks(msRaw, scan)
              if(!is.na(deprofile))
              {								
                pks <- deprofile.scan(
                    pks, method = deprofile, noise = NA, colnames = FALSE
                )
              }
              colnames(pks) <- c("mz","int")
              return(pks)
            }
        )
        return(list(
                foundOK = TRUE,
                parentScan = ret$parentScan,
                parentHeader = as.data.frame(ret$parentHeader),
                childScans = childScans,
                childHeaders= childHeaders,
                parentPeak=ret$parentPeak,
                peaks=msmsPeaks
            #xset=xset#,
            #msRaw=msRaw
            ))
      })
  names(spectra) <- eic$acquisitionNum
  return(spectra)
}
# 
# Author: stravsmi
###############################################################################
msmsRead.ticms2 <- function(w, filetable = NULL, files = NULL, cpdids = NULL, 
                     readMethod, mode, confirmMode = FALSE, useRtLimit = TRUE, 
                     Args = NULL, settings = getOption("RMassBank"), progressbar = "progressBarHook", MSe = FALSE){
  
  ##Read the files and cpdids according to the definition
  ##All cases are silently accepted, as long as they can be handled according to one definition
  if(is.null(filetable)){
    ##If no filetable is supplied, filenames must be named explicitly
    if(is.null(files))
      stop("Please supply the files")
    
    ##Assign the filenames to the workspace
    w@files <- unlist(files)
    
    ##If no filetable is supplied, cpdids must be delivered explicitly or implicitly within the filenames
    if(is.null(cpdids)){
      splitfn <- strsplit(files,"_")
      splitsfn <- sapply(splitfn, function(x) x[length(x)-1])
      if(suppressWarnings(any(is.na(as.numeric(splitsfn)[1]))))
        stop("Please supply the cpdids corresponding to the files in the filetable or the filenames")
      cpdids <- splitsfn
    }
  } else{
    ##If a filetable is supplied read it
    tab <- read.csv(filetable, stringsAsFactors = FALSE)
    w@files <- tab[,"Files"]
    cpdids <- tab[,"ID"]
  }
  
  ##If there's more cpdids than filenames or the other way around, then abort
  if(length(w@files) != length(cpdids)){
    stop("There are a different number of cpdids than files")
  }
  if(!(readMethod %in% c("mzR","peaklist","xcms","ticms2"))){
    stop("The supplied method does not exist")
  }
  if(!all(file.exists(w@files))){
    stop("The supplied files don't exist")
  }
  
  ##This should work
  if(readMethod == "mzR"){
    ##Progressbar
    nLen <- length(w@files)
    nProg <- 0
    pb <- do.call(progressbar, list(object=NULL, value=0, min=0, max=nLen))
    
    count <- 1
    w@specs <-  lapply(w@files, function(fileName){
      spec <- findMsMsHR(fileName, cpdids[count], mode, confirmMode, useRtLimit,
                         ppmFine = settings$findMsMsRawSettings$ppmFine,
                         mzCoarse = settings$findMsMsRawSettings$mzCoarse,
                         fillPrecursorScan = settings$findMsMsRawSettings$fillPrecursorScan,
                         rtMargin = settings$rtMargin,
                         deprofile = settings$deprofile)
      
      ## Progress:
      nProg <<- nProg + 1
      pb <- do.call(progressbar, list(object=pb, value= nProg))
      
      ##Counting the index of cpdids
      count <<- count + 1
      return(spec)
    })
    names(w@specs) <- basename(as.character(w@files))
    return(w)
  }
  
  ##This should work
  if(readMethod == "ticms2"){
    ##Progressbar
    nLen <- length(w@files)
    nProg <- 0
    pb <- do.call(progressbar, list(object=NULL, value=0, min=0, max=nLen))
    
    count <- 1
    w@specs <-  lapply(w@files, function(fileName){
      spec <- findMsMsHR.ticms2.file(fileName, cpdids[count], mode, confirmMode, useRtLimit,
                                     ppmFine = settings$findMsMsRawSettings$ppmFine,
                                     mzCoarse = settings$findMsMsRawSettings$mzCoarse,
                                     fillPrecursorScan = settings$findMsMsRawSettings$fillPrecursorScan,
                                     rtMargin = settings$rtMargin,
                                     deprofile = settings$deprofile)
      
      ## Progress:
      nProg <<- nProg + 1
      pb <- do.call(progressbar, list(object=pb, value= nProg))
      
      ##Counting the index of cpdids
      count <<- count + 1
      return(spec)
    })
    names(w@specs) <- basename(as.character(w@files))
    return(w)
  }
  
  ##Peaklist-readmethod 
  if(readMethod == "peaklist"){
    w <- createSpecsFromPeaklists(w, cpdids, filenames=w@files, mode=mode)
    names(w@specs) <- basename(as.character(w@files))
    return(w)
  }
  
  ##xcms-readmethod 
  if(readMethod == "xcms"){
    ufiles <- unique(w@files)
    uIDs <- unique(cpdids)
    ##Routine for the case of multiple cpdIDs per file and multiple files per cpdID
    dummySpecs <- list()
    w@specs <- list()
    for(i in 1:length(ufiles)){ ##Create list
      dummySpecs[[i]] <- newMsmsWorkspace()
      dummySpecs[[i]]@specs <- list()
      FileIDs <- cpdids[which(w@files == ufiles[i])]
      metaSpec <- findMsMsHRperxcms.direct(ufiles[i], FileIDs, mode=mode, findPeaksArgs=Args, MSe = MSe)
      for(j in 1:length(FileIDs)){
        dummySpecs[[i]]@specs[[length(dummySpecs[[i]]@specs)+1]] <- metaSpec[[j]]
      }
      
    }
    
    if(length(dummySpecs) > 1){
      for(j in 2:length(dummySpecs)){
        dummySpecs[[1]] <- c.msmsWSspecs(dummySpecs[[1]],dummySpecs[[j]])
      }
    }
    
    ##You need as many names as there were different IDs
    ##And the Names and IDs have to go together in some way
    ##Find out Names that make sense: (cpdID with Name of File that uses cpdID)
    FNames <- vector()
    for(i in uIDs){
      nindex <- min(which(i == cpdids))
      FNames <- c(FNames,paste(w@files[nindex],"_",cpdids[nindex],sep=""))
    }
    
    w@specs <- dummySpecs[[1]]@specs
    names(w@specs) <- basename(as.character(FNames))
    w@files <- basename(as.character(FNames))
    return(w)
  }
}

findMsMsHR.ticms2.d <- function(msRaw, cpdID, mode = "pH", confirmMode = 0, useRtLimit = TRUE, 
    ppmFine = getOption("RMassBank")$findMsMsRawSettings$ppmFine,
    mzCoarse = getOption("RMassBank")$findMsMsRawSettings$mzCoarse,
    fillPrecursorScan = getOption("RMassBank")$findMsMsRawSettings$fillPrecursorScan,
    rtMargin = getOption("RMassBank")$rtMargin,
    deprofile = getOption("RMassBank")$deprofile,
    headerCache = NA)
{
  # for finding the peak RT: use the gauss-fitted centwave peak
  # (centroid data converted with TOPP is necessary. save as
  # mzData, since this is correctly read :P)
  #xset <- xcmsSet(fileName, method="centWave",ppm=5, fitgauss=TRUE)
  
  # find cpd m/z
  mzLimits <- findMz(cpdID, mode)
  mz <- mzLimits$mzCenter
  limit.fine <- ppm(mz, ppmFine, p=TRUE)
  if(!useRtLimit)
    rtLimits <- NA
  else
  {
    dbRt <- findRt(cpdID)
    rtLimits <- c(dbRt$RT - rtMargin, dbRt$RT + rtMargin) * 60
  }
  spectra <- findMsMsHR.ticms2(msRaw, mz, mzCoarse, limit.fine, rtLimits, confirmMode + 1,headerCache
      ,fillPrecursorScan, deprofile)
  # check whether a) spectrum was found and b) enough spectra were found
  if(length(spectra) < (confirmMode + 1))
    sp <- list(foundOK = FALSE)
  else
    sp <- spectra[[confirmMode + 1]]
  
  sp$mz <- mzLimits
  sp$id <- cpdID
  sp$formula <- findFormula(cpdID)
  return(sp)
}

findMsMsHR.ticms2.file <- function(fileName, cpdID, mode="pH",confirmMode =0, useRtLimit = TRUE,
    ppmFine = getOption("RMassBank")$findMsMsRawSettings$ppmFine,
    mzCoarse = getOption("RMassBank")$findMsMsRawSettings$mzCoarse,
    fillPrecursorScan = getOption("RMassBank")$findMsMsRawSettings$fillPrecursorScan,
    rtMargin = getOption("RMassBank")$rtMargin,
    deprofile = getOption("RMassBank")$deprofile)
{
  
  # access data directly for finding the MS/MS data. This is done using
  # mzR.
  msRaw <- openMSfile(fileName)
  ret <- findMsMsHR.ticms2.d(msRaw, cpdID, mode, confirmMode, useRtLimit, ppmFine, mzCoarse, fillPrecursorScan,
      rtMargin, deprofile)
  mzR::close(msRaw)
  return(ret)
}
