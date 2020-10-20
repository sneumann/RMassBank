[![Build Status](https://travis-ci.org/MassBank/RMassBank.svg?branch=main)](https://travis-ci.org/MassBank/RMassBank)

# RMassBank

Workflow to process tandem MS files and build MassBank records. Functions include automated extraction of tandem MS spectra, formula assignment to tandem MS fragments, recalibration of tandem MS spectra with assigned fragments, spectrum cleanup, automated retrieval of compound information from Internet databases, and export to MassBank records.

Authors: Michael Stravs, Emma Schymanski, Steffen Neumann, Erik Mueller, with contributions from Tobias Schulze

Maintainer: RMassBank at Eawag <massbank at eawag.ch>

Citation (from within R, enter `citation("RMassBank")`):

Stravs MA, Schymanski EL, Singer H and Hollender J (2013). “Automatic Recalibration and Processing of Tandem Mass Spectra using Formula Annotation.” Journal of Mass Spectrometry, 48(1), pp. 188.


# Branch and merge policy

We aim to have a `main` branch that is in sync with BioC `master` and always passes the Travis CI checks.
All development should take place in the `dev` branch and via Pull Requests.

Note: to push towards BioC you can `git checkout master` (which is the BioC `master`), then merge the github branch via `git merge main` and `git push upstream master` (assuming the BioC remote is called `upstream` as recommended).
