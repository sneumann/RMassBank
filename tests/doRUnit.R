#### doRUnit.R --- Run RUnit tests
####------------------------------------------------------------------------

### Structure borrowed from rcppgls:
### https://github.com/eddelbuettel/rcppgsl/blob/master/tests/doRUnit.R

if(require("RUnit", quietly = TRUE)) {
	if(require("RMassBankData", quietly = TRUE)) {
		pkg <- "RMassBank"
		print("Starting tests")
		require(pkg, character.only=TRUE)

		path <- system.file("unitTests", package = pkg)

		stopifnot(file.exists(path), file.info(path.expand(path))$isdir)

		source(file.path(path, "runTests.R"), echo = TRUE)
	} else {
		print("Package RMassBankData not available, cannot run unit tests")
	}
} else {
	print("Package RUnit not available, cannot run unit tests")
}            