
#' Minimalistic LOCF and LOCB functions.
#' 
#' Thanks to http://stackoverflow.com/a/2776751 for the idea.
#' But this here is recoded using which() and therefore not CC-BY-SA :)
#' (Also, it fills the first unmatcheds with NA and not with the first value)
#' Yes, I know that library(zoo) has na.locf, but zoo is GPL-2.
#' 
#' @aliases .locf .locb
#' 
#' @noRd 
#' @examples
#' vec <- c(NA, NA,  "A", "B", NA, NA, NA,  "C", NA,  "D", NA, NA, NA)
#' .locb(vec)
#' .locf(vec)
#' 
.locf <- function(vec) {
	vec.isna <- !is.na(vec)
	vec.isna.index <- cumsum(vec.isna)
	# The construct below deals with the fact that vec[0] == integer(0) and not NA...
	# Therefore we retrieve indexes to c(NA,vec) from c(0, which(vec.isna)) + 1
	# (the +1 shifts all values so the 0s point to NA).
	c(NA,vec)[c(0,which(vec.isna))[vec.isna.index+1]+1]
}
.locb <- function(vec) {
	rvec <- rev(vec)
	rev(.locf(rvec))
}


#' De-profile a high-resolution MS scan in profile mode.
#' 
#' The \code{deprofile} functions convert profile-mode high-resolution input data to "centroid"-mode
#' data amenable to i.e. centWave. This is done using full-width, half-height algorithm, spline
#' algorithm or local maximum method.
#' 
#' The \code{deprofile.fwhm} method is basically an R-semantic version of the "Exact Mass" m/z deprofiler
#' from MZmine. It takes the center between the m/z values at half-maximum intensity for the exact peak mass.
#' The logic is stolen verbatim from the Java MZmine algorithm, but it has been rewritten to use the fast
#' R vector operations instead of loops wherever possible. It's slower than the Java implementation, but 
#' still decently fast IMO. Using matrices instead of the data frame would be more memory-efficient
#' and also faster, probably.
#' 
#' The \code{deprofile.localMax} method uses local maxima and is probably the same used by e.g. Xcalibur.
#' For well-formed peaks, "deprofile.fwhm" gives more accurate mass results; for some applications,
#' \code{deprofile.localMax} might be better (e.g. for fine isotopic structure peaks which are
#' not separated by a valley and also not at half maximum.) For MS2 peaks, which have no isotopes,
#' \code{deprofile.fwhm} is probably the better choice generally.
#' 
#' \code{deprofile.spline} calculates the mass using a cubic spline, as the HiRes peak detection
#' in OpenMS does.
#' 
#' The word "centroid" is used for convenience to denote not-profile-mode data.
#' The data points are NOT mathematical centroids; we would like to have a better word for it.
#' 
#' The \code{noise} parameter was only included for completeness, I personally don't use it.
#' 
#' \code{deprofile.fwhm} and \code{deprofile.localMax} are the workhorses; 
#' \code{deprofile.scan} takes a 2-column scan as input.
#' \code{deprofile} dispatches the call to the appropriate worker method.
#' 
#' @note Known limitations: If the absolute leftmost stick or the absolute rightmost stick in
#' 				a scan are maxima, they will be discarded! However, I don't think this will
#' 				ever present a practical problem.
#' 
#' @aliases deprofile.scan deprofile.fwhm deprofile.localMax deprofile.spline
#' @usage		deprofile.scan(scan, noise = NA, method = "deprofile.fwhm", 
#' 					colnames = TRUE, ...)
#' 
#' 				deprofile(df, noise, method, ...)
#' 
#' 				deprofile.fwhm(df, noise = NA, cut = 0.5)
#' 
#' 				deprofile.localMax(df, noise = NA)
#' 
#' 				deprofile.spline(df, noise=NA, minPts = 5, step = 0.00001)
#' @param scan A matrix with 2 columns for m/z and intensity; from profile-mode high-resolution data. Corresponds
#' 				to the spectra obtained with xcms::getScan or mzR::peaks. 
#' @param noise The noise cutoff. A peak is not included if the maximum stick intensity of the peak
#' 				is below the noise cutoff.
#' @param method "deprofile.fwhm" for full-width half-maximum or "deprofile.localMax" for
#' 				local maximum.
#' @param colnames For deprofile.scan: return matrix with column names (xcms-style, 
#' 				\code{TRUE}, default) or plain (mzR-style, \code{FALSE}).
#' @param df A dataframe with at least the columns \code{mz} and \code{int} to perform deprofiling on. 
#' @param ... Arguments to the workhorse functions \code{deprofile.fwhm} etc.
#' @param cut A parameter for \code{deprofile.fwhm} indicating where the peak flanks should be taken. Standard is 0.5
#' 				(as the algorithm name says, at half maximum.) Setting \code{cut = 0.75} would instead determine the peak
#' 				width at 3/4 maximum, which might give a better accuracy for merged peaks, but could be less accurate
#' 				if too few data points are present.
#' @param minPts The minimal points count in a peak to build a spline; for peaks with less
#' 				points the local maximum will be used instead. Note: The minimum value
#' 				is 4!
#' @param step The interpolation step for the calculated spline, which limits the maximum 
#' 				precision which can be achieved. 
#' @return \code{deprofile.scan}: a matrix with 2 columns for m/z and intensity
#' 
#' @examples 
#' \dontrun{
#' mzrFile <- openMSfile("myfile.mzML")
#' acqNo <- xraw@@acquisitionNum[[50]]
#' scan.mzML.profile <- mzR::peaks(mzrFile, acqNo)
#' scan.mzML <- deprofile.scan(scan.mzML.profile) 
#' close(mzrFile)
#' }
#' 
#' @author Michael Stravs, Eawag <michael.stravs@@eawag.ch>
#' @references 
#' mzMine source code \href{http://sourceforge.net/svn/?group_id=139835}{http://sourceforge.net/svn/?group_id=139835}
#' @export
deprofile <- function(df, noise, method, ...)
{
	return(do.call(method, list(df, noise, ...)))
}

#' @export
deprofile.scan <- function(scan, noise = NA, method="deprofile.fwhm", colnames = TRUE, ...)
{
	# Format the input
	df <- as.data.frame(scan)
	colnames(df) <- c("mz", "int")
	# Call the actual workhorse
	peaklist <- deprofile(df, noise, method, ...)
	# return the centroided peaklist
	peaklist.m <- as.matrix(peaklist[,c("mz", "int")])
	if(colnames)
		colnames(peaklist.m) <- c("mz", "intensity")
	else
		colnames(peaklist.m) <- NULL
	return(peaklist.m)
}

#' @export
deprofile.fwhm <- function(df, noise=NA, cut=0.5)
{
	# split sticks into groups according to how MzMine does it:
	# a new group starts at zeroes and at new ascending points
	df$groupmax <- NA
	rows <- nrow(df)
	df <- within(df,
			{
				# identify local maxima
				groupmax[which(diff(sign(diff(df$int)))<0) + 1] <- which(diff(sign(diff(df$int)))<0) + 1
				# make forward-filled and backward-filled list for which was the last maximum.
				# This assigns the sticks to groups.
				groupmax_f <- .locf(groupmax)
				groupmax_b <- .locb(groupmax)
				# take backward-filled group as default
				# and forward-filled group where the peak was ascending
				groupmax_b[which(diff(df$int)<0)+1 ] <- groupmax_f[which(diff(df$int) <0)+1]
				# eliminate zeroes
				groupmax_b[which(df$int==0)]<-NA
				# add "next intensities" and "next m/z" as well as "index" (n)
				n <- 1:rows
				int1 <- c(df$int[-1],0)
				mz1 <- c(df$mz[-1],0)
				# delete all the intensity+1 values from points which are last-in-group and therefore have int1 from next group!
				int1[which(!is.na(groupmax))-1] <- 0
				# find maximal intensity point for each peak member
				maxint <- df[groupmax_b, "int"]
				hm <- maxint * cut
				up <- ifelse( (df$int <= hm) & (int1 >= hm) & (n < groupmax_b), groupmax_b, NA)
				down <- ifelse( (df$int >= hm) & (int1 <= hm) & (n > groupmax_b), groupmax_b, NA)
			})
	# Compile finished peak list
	peaklist <- df[which(!is.na(df$groupmax)),]
	# Noise parameter:
	if(!is.na(noise))
		peaklist <- peaklist[peaklist$maxint > noise,]
	# Find which peaks might have a FWHM value to substitue for the maxint mz value
	# We isolate the peaklists for left-hand and for right-hand FWHM peak.
	# If any straight line is found, the rightmost of the left points and the leftmost of the right points is used.
	peaklist.left <- df[which(!is.na(df$up) & !duplicated(df$up, fromLast=TRUE)),]
	peaklist.right <- df[which(!is.na(df$down) & !duplicated(df$down)),]
	# calculate the slopes and the corresponding m/z value at half maximum
	peaklist.left$slope <- (peaklist.left$int1 - peaklist.left$int) / (peaklist.left$mz1 - peaklist.left$mz)
	peaklist.right$slope <- (peaklist.right$int1 - peaklist.right$int) / (peaklist.right$mz1 - peaklist.right$mz)
	peaklist.left$mzleft <- peaklist.left$mz + (peaklist.left$hm - peaklist.left$int) / peaklist.left$slope
	peaklist.right$mzright <- peaklist.right$mz + (peaklist.right$hm - peaklist.right$int) / peaklist.right$slope
	# add the two values to the full-peaklist where they exist
	peaklist <- merge(peaklist, peaklist.left[,c("up", "mzleft")], by.x="groupmax", by.y="up", all.x=TRUE, suffix=c("", ".left"))
	peaklist <- merge(peaklist, peaklist.right[,c("down", "mzright")], by.x="groupmax", by.y="down", all.x=TRUE, suffix=c("", ".right"))
	# Find which entries have both a left and a right end,
	# and calculate the center mass for them.
	peaklist.indexMzhm <- which(!is.na(peaklist$mzleft) & !is.na(peaklist$mzright))
	peaklist[peaklist.indexMzhm, "mz"] <- (peaklist[peaklist.indexMzhm, "mzleft"] + peaklist[peaklist.indexMzhm, "mzright"]) / 2
	
	return(peaklist)
}

#' @export
deprofile.localMax <- function(df, noise=NA)
{
	# split sticks into groups:
	# a new group starts at zeroes and at new ascending points
	df$groupmax <- NA
	rows <- nrow(df)
	df$groupmax[which(diff(sign(diff(df$int)))<0) + 1] <- which(diff(sign(diff(df$int)))<0) + 1
	peaklist <- df[which(!is.na(df$groupmax)),]
	# Noise parameter:
	if(!is.na(noise))
		peaklist <- peaklist[peaklist$int > noise,]
	# And that's it.
	return(peaklist)
}

# This spline thing will be very slow :)
#' @export
deprofile.spline <- function(df, noise=NA, minPts = 5, step= 0.00001)
{
	df$groupmax <- NA
	rows <- nrow(df)
	# Group the peaks like the FWHM routine
	df <- within(df,
			{
				# identify local maxima
				groupmax[which(diff(sign(diff(df$int)))<0) + 1] <- which(diff(sign(diff(df$int)))<0) + 1
				# make forward-filled and backward-filled list for which was the last maximum.
				# This assigns the sticks to groups.
				groupmax_f <- .locf(groupmax)
				groupmax_b <- .locb(groupmax)
				# take backward-filled group as default
				# and forward-filled group where the peak was ascending
				groupmax_b[which(diff(df$int)<0)+1 ] <- groupmax_f[which(diff(df$int) <0)+1]
				# eliminate zeroes
				groupmax_b[which(df$int==0)]<-NA
			})
	groups <- na.omit(df$groupmax)
	peaklist <- t(sapply(groups, function(group)
			{
				pk <- df[which(df$groupmax_b == group),]
				# if there are not enough points, return the local maximum
				if(nrow(pk) < minPts)
					return(as.matrix(pk[which(pk$groupmax == group),c("mz", "int")]))
				# fit a spline
				spl <- smooth.spline(pk$mz, pk$int)
				# predict in small steps
				pred <- seq(from=min(pk$mz), to=max(pk$mz), by=step)
				curve <- predict(spl, x=pred)
				# find top and return it as peak
				top <- which.max(curve$y)
				return(as.matrix(c(curve$x[[top]], curve$y[[top]])))
			}))
	colnames(peaklist) <- c("mz", "int")
	peaklist[,"mz"] <- unlist(peaklist[,"mz"])
	peaklist[,"int"] <- unlist(peaklist[,"int"])
	return(peaklist)
}
