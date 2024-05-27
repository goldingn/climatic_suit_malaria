# Load libraries
library(sp)
library(raster)

# Load a single layer as an example
temp_stack <- stack("AllData/2022_03.grib")
sample_layer <- temp_stack[[1]]

# Calculate which coordinates are not NA 
basena = extract(calc(sample_layer, function(x)(is.na(x))), 1:ncell(sample_layer))
valid_cells = which(basena!=1)

# Generate blank grd for plotting later
blank <- temp_stack[[1]]
blank[valid_cells] <- 0
writeRaster(blank, filename="blankNew.grd",overwrite=TRUE)

# We run a degree day accumulation across the summer months
degree_day_accumulation = function() {
  file_names <- c("AllData/2022_06.grib",
                  "AllData/2022_07.grib",
                  "AllData/2022_08.grib",
                  "AllData/2022_12.grib",
                  "AllData/2023_01.grib",
                  "AllData/2023_02.grib")
  current_layer = 1
  in_place_accumulation = numeric(length(valid_cells))
  for(fileName in file_names) {
    temp_brick = brick(fileName)
    number_of_layers = nlayers(temp_brick)
    layers = seq(2, number_of_layers - 2, 4)
    for (layer in layers) {
      brick_layer = pmax((temp_brick[[layer]][valid_cells] - 289.15) / 12, 0)
      in_place_accumulation = in_place_accumulation + brick_layer
      print(paste(fileName, toString(layer), "Current Layer", current_layer, sep=" "))
      current_layer = current_layer + 1
    }
  }
  return(in_place_accumulation)
}
accumulation = degree_day_accumulation()
# Cells need 111 degree days to have a chance at producing a viable population
# Since this accumulation covers half the runtime of the model, we consider only cells with 55 or more degree days
# These cells will be excluded from the model to reduce its running time
accumulated_indices = which(accumulation >= 55)
valid_cells = valid_cells[accumulated_indices]

# Split valid cells into chunks for running
number_of_chunks = 25
number_of_cells = length(valid_cells)
cells_per_chunk = ceiling(number_of_cells / number_of_chunks)
valid_cells_indices = seq(1, number_of_cells, cells_per_chunk)
valid_cells_matrix = matrix(NA, number_of_chunks, cells_per_chunk)

# For loop fills in all rows except the last.
# Last row is filled in manually afterwards as it has a different number of elements
for(chunk in 1:(number_of_chunks - 1)) {
  valid_cells_matrix[chunk,] = valid_cells[valid_cells_indices[chunk]:(valid_cells_indices[chunk + 1] - 1)]
}
cells_in_final_chunk = number_of_cells - valid_cells_indices[number_of_chunks] + 1
valid_cells_matrix[number_of_chunks,1:cells_in_final_chunk] = 
  valid_cells[valid_cells_indices[number_of_chunks]:number_of_cells]

save(valid_cells_matrix, file="ValidCells.RData")
