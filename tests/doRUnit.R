#### doRUnit.R --- Run RUnit tests
####------------------------------------------------------------------------

### Structure borrowed from rcppgls:
### https://github.com/eddelbuettel/rcppgsl/blob/master/tests/doRUnit.R

if(require("RUnit", quietly = TRUE)) {
	if(require("RMassBankData", quietly = TRUE) && !(compareVersion(installed.packages()["RMassBankData","Version"],"1.21.0") == -1)) {
		pkg <- "RMassBank"
		print("Starting tests")
		require(pkg, character.only=TRUE)

		path <- system.file("unitTests", package = pkg)

		stopifnot(file.exists(path), file.info(path.expand(path))$isdir)

		source(file.path(path, "runTests.R"), echo = TRUE)
	} else {
		## Taking this message out until the new RMassBankData is on bioc, just to avoid confusion.
        # message("Package RMassBankData with version > 1.99 not available, cannot run unit tests")
	}
} else {
	message("Package RUnit not available, cannot run unit tests")
}
