In process of refactoring to work with ERA5 data. Currently working with .grib files downloaded directly from https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land?tab=form
Each grib file covers one month in 2 hr increments (00:00, 02:00,...22:00), and is about 2GB, so haven't uploaded to the repo for now.

Running with the expanded dataset will be a challenge - R cannot handle all 5000 raster layers at once, so each loop of the model now includes the data extraction step.

* TempModel.slurm - this should run all 25 chunks of the Temperature model
* TempModel_Run.R - this extracts data to a file TempMatrixN.RData (where N is the chunk run), and runs the model
* TempModel_Initial.R - this file generates the necessary ValidCells.RData file, can be safely ignored assuming this RData file is present and up to date
* ValidCells.RData - the temp model will attempt to load this file as part of its run

* HumModel.slurm - this should run all 25 chunks of the Humidity/Temperature model
* HumModel_Run.R - this extracts data to a file TempDewpointMatrixN.RData, and runs the model. Requires the file TempMatrixN.RData to be loaded

* RainfallModel.slurm - this should run all 25 chunks of the Rainfall/Humidity/Temperature model
* RainfallModel_Run.R - this extracts data to a file FullMatrixN.RData, and runs the model. Requires the file TempDewpointMatrixN.RData to be loaded

