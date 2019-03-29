# gaiadr2-6dgold-example
Example of the construction of a sample of sources from Gaia DR2 with all recommend data quality filtering applied.

## Goals of this example project
* Show how to contruct from the Gaia DR2 archive data a 'clean', aka 'gold', sample of sources for 6D phase
  space studies to which all filtering recommendations from the Gaia Collaboration have been applied.
* The resulting data set is intended for those who want to start working on Gaia DR2 data without having
  to worry about data quality issues (such as spurious astrometry, radial velocities, unreliable
  photometry, etc). 
  * Undergraduate students who want to do a Gaia research project.
  * The gold sample as a teaching tool; what should I expect to see in a clean Gaia DR2 sample?
* The 'gold' sample is thus 'pure' but certainly not complete!

## Installation
* Clone the repository or download the Python notebook.
* You must download the following files and (after unzipping if relevant) and drop them in the `data` folder:
  * [Tables with RUWE normalization factors](https://www.cosmos.esa.int/documents/29201/1769576/DR2_RUWE_V1.zip/d90f37a8-37c9-81ba-bf59-dd29d9b1438f)
  * [List of Gaia DR2 sources with potentially spurious radial velocities](https://arxiv.org/src/1901.10460v1/anc/rvscontamination.csv)
  * Download the necessary input data from the Gaia DR2 archive using the ADQL queries in the notebook
    and drop the result in the  `data` folder. __NOTE__ FITS format is used in this example, so make sure
    to request the download tables as FITS.

## Required python packages
* [Numpy](https://www.numpy.org/)
* [Matplotlib](https://matplotlib.org/)
* [Astropy](http://www.astropy.org/)
* [ruwetools](https://github.com/agabrown/gaiadr2-ruwe-tools) (Already included for convenience.)

## TODO
* Document limitations of the gold sample
  * Examples of the selection effects introduced when filtering on data quality indicators.
  * Warning on incorrect use of parallaxes.
* Error propagation for the Galactocentric phase space coordinates.
