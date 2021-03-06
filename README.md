# gaiadr2-6dgold-example
Example of the construction of a sample of sources from Gaia DR2 with all recommend data quality filtering applied.

## Credits
This project grew out of the [Santa Barbara Gaia Sprint 2019](http://gaia.lol/2019SB.html), hosted at KITP, and the KITP programme ["Dynamical Models for Stars and Gas in Galaxies in the Gaia Era"](https://www.kitp.ucsb.edu/activities/gaia19). Thanks to Juna Kollmeier for suggesting the creation of Gaia Gold samples.

## Goals of this example project
* Show how to contruct from the Gaia DR2 archive data a 'clean', aka 'gold', sample of sources for 6D phase
  space studies to which all data quality control recommendations from the Gaia Collaboration have been applied.
  * In particular the recommendations from the [Gaia DR2 Known Issues pages](https://www.cosmos.esa.int/web/gaia/dr2-known-issues).
* The resulting data set is intended for those who want to start working on Gaia DR2 data without having
  to worry about data quality issues (such as spurious astrometry, radial velocities, unreliable
  photometry, etc).
  * Undergraduate students who want to do a Gaia research project.
  * The gold sample as a teaching tool; what should I expect to see in a clean Gaia DR2 sample?
* The 'gold' sample is thus 'pure' but certainly not complete and will suffer from non-trivial selection
  effects!

## Installation
* Clone the repository or download the Python notebook.
* You must download the following files and (after unzipping if relevant) and drop them in the `data` folder:
  * [Tables with RUWE normalization factors](https://www.cosmos.esa.int/documents/29201/1769576/DR2_RUWE_V1.zip/d90f37a8-37c9-81ba-bf59-dd29d9b1438f)
  * [List of Gaia DR2 sources with potentially spurious radial velocities](https://arxiv.org/src/1901.10460v1/anc/rvscontamination.csv)
  * Download the necessary input data from the Gaia DR2 archive using the ADQL queries in the notebook
    and drop the result in the  `data` folder. __NOTE__ FITS format is used in this example, so make sure
    to request the download tables as FITS.

## Usage
* The notebook explains how the 6D Gold sample is created from the data downloaded from the Gaia DR2
  archive.
* Use the command-line tool `build_6dgold_tables.py` to build the 6D Gold table(s) in FITS format so the
  Gold sample be re-used easily, loaded into [topcat](http://www.star.bris.ac.uk/~mbt/topcat/), etc.
* __NOTE__ that in the notebook as well as in the command line tool specific choices are made for the Milky
  parameters (such as the distance from the Sun to the Galactic centre) in the calculation of the Galactocentric
  quantities. These choices can be adapted in the notebook and in the `ggtools.py` code.

## Required python packages
* [Numpy](https://www.numpy.org/)
* [Matplotlib](https://matplotlib.org/)
* [Astropy](http://www.astropy.org/)
* [ruwetools](https://github.com/agabrown/gaiadr2-ruwe-tools) (Already included for convenience.)

## TODO
* Make the resulting FITS tables available somewhere.
* Document limitations of the gold sample
  * Examples of the selection effects introduced when filtering on data quality indicators.
  * Warning on incorrect use of parallaxes.
* Error propagation for the Galactocentric phase space coordinates.
