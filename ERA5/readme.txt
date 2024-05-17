In process of refactoring to work with ERA5 data. Currently working with .grib files downloaded directly from https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land?tab=form
Each grib file covers one month in 2 hr increments (00:00, 02:00,...22:00), and is about 2GB, so haven't uploaded to the repo for now.

Running with the expanded dataset will be a challenge - R cannot handle all 5000 raster layers at once, so each loop of the model now includes the data extraction step.

Files:
* TempModel.slurm - this should run all 25 chunks of the model
* TempModel_Run.R - this file contains the model as well as the data extraction logic
* TempModel_Run_Fast.R - this file just contains the model, assumes data has been extracted and saved already
* TempModel_Initial.R - this file generates the necessary ValidCells.RData file, can be safely ignored assuming this RData file is present and up to date
* ValidCells.RData - the temp model will attempt to load this file as part of its run
