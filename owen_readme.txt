Keeping this here as a reference on how to use this collection.

* Existing MonthlyData files can be run through a block in Rasterize.R, to "Take monthly data and convert into Rasters". This will give 12 .gri files, one for each month.
* If running on temp model, import ProjectNew.RData. Otherwise, you will want to import HProjectNew.RData

* The model can also be re-run. Slurm scripts are provided to run the program on the High Performance computing system at Unimelb.
* To re-run, the script (H|S)ProjectInputNew.R and dependent data file SInputNew.RData are needed. 
* The R script will change depending on which model you want to be run (S prefix denotes with rainfall, humidity and temp, H prefix denotes only with humidity and temp, no prefix is with only temp)

* The script (H|S)ProjectBuildNew is what generates the dependent data file. The order 'no prefix' -> H -> S gradually builds up the RData file, with each script adding another block of data
* These scripts will be dependent on temperature, humidity, rainfall and altitude data, which might need to be standardised to fit with the 5km resolution the scripts are designed for