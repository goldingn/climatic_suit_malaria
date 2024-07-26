# try downloading era5 data and plugging into Owen's analysis code

# Note: there is an mcera5 package to access data from the copernicus portal,
# but its dependencies on microclima and the deprecated rgdal package make
# installation difficult.
# It is also possible to access ERA5 from google earth engine using the rgee
# package, but registration and access (via python) is similarly painful.

# Therefore I am instead making a manual data request in NetCDF format via
# copernicus.eu, and setting up code to process the downloaded file, using
# code in the mcera package.

# download from https://cds.climate.copernicus.eu, after registering.

# download the following:

# Product type: Monthly averaged reanalysis by hour of day
# Variables:
#     2m dewpoint temperature,
#     2m temperature,
#     Surface pressure,
#     Total precipitation
# Years: 2014-2024 (inclusive)
# Months: All
# Times: All (hourly 00:00-23:00)
# Region: Whole available region:
# Format: NetCDF (experimental)

# These variables can be used to obtain precipitation, temperature, and humidity
# humidity can be calculated from dewpoint temperature, air temperature, and
# pressure using:
mcera5:::humfromdew()

# The ERA5 netcdf values can be subset to specific lat-longs (and times) and
# converted to more useful units using mcera5::extract_clim() (which calls
# mcera5:::nc_to_df()) or manually with tidync & tidyverse functions below.

# nc_filename <- "?"
# # load the netcdf file
# nc_filename %>%
#   ncdf4::nc_open() %>%
#   # connect to tidync
#   tidync::tidync() %>%
#   # prepare filters and to be used in the hyper_tibble call
#   # tidync::hyper_filter(...) %>%
#   # apply the filter to read in the data as a tibble
#   tidync::hyper_tibble() %>%
#   # do manipulation on this object 
  
