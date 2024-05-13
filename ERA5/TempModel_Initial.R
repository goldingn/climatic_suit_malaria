# Load libraries
library(sp)
library(raster)

# Load a single layer as an example
temp_stack <- stack("TempData/2m_temp_2022_03.grib")
sample_layer <- temp_stack[[1]]

# Calculate which coordinates are not NA 
basena = extract(calc(sample_layer, function(x)(is.na(x))), 1:ncell(sample_layer))
valid_cells = which(basena!=1)

# Generate blank grd for plotting later
blank <- temp_stack[[1]]
blank[valid_cells] <- 0
writeRaster(blank, filename="blankNew.grd",overwrite=TRUE)

# TODO: Do we want to do further filtering to cells which never get a high temperature?

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
