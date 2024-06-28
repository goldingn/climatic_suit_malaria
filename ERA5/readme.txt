In order to run this code, you will need the raw grib fies downloaded directly from https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land?tab=form.
Request a file with 2m dewpoint temperature, 2m temperature, Evaporation from open water surfaces excluding oceans and Total precipitation.
All days should be selected, and 2hr increments of time should be selected (00:00, 02:00,...22:00). For our project, the whole available region was selected, although the code should work with sub-regions.
Once downloaded, rename the file to be "yyyy_mm.grib", using the year and month of the that particular file, and placed in a folder called 'AllData'
This code is hard coded to take the months from March 2022 to May 2023. The model can be easily refactored to take a different date range, although some of the utils functions for plotting results may need rewriting.

* TempModel_Initial.R - this file generates the necessary ValidCells.RData file, and the blank raster for plotting. These files are included in the repo, but would need to be regenerated if running the model on a different region.
* TempModel_Run.R - this extracts data and runs the temperature-only model
* HumModel_Run.R - this extracts data and runs the temperature-humidity model. In order to speed up the model, this file assumes the data extraction in TempModel has already been done.
* RainfallModel_Run.R - this extracts data and runs the temperature-humidity-rainfall model. In order to speed up the model, this file assumes the data extraction in HumModel has already been done.
* Utils.R - this file contains utility functions for mapping model results to output images

The slurm files can be used to the run the models through sbatch.
For our model we split the workload into 25 chunks for better parallel processing - this could be done with more or fewer, but would require the Initial file to be re-run, and some utils may need refactoring.

