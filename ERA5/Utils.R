convertCellNumberIntoCoords = function(cell) {
  offsetCell = cell - 1
  xcoord = 90 - (floor(offsetCell / 3600) / 10)
  ycoord = ((offsetCell %% 3600) / 10) - 180
  return (c(xcoord, ycoord))
}

convertCoordsIntoCellNumber = function(coords) {
  numberCellsDown = (90 - coords[1]) * 10
  numberCellsAcross = (coords[2] + 180) * 10
  return ((numberCellsDown * 3600) + numberCellsAcross)
}

