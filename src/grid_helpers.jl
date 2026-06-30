_spatial_grid(grid) = grid
_spatial_grid(spectral_grid::SpectralGrid) = spectral_grid.grid

_grid_npoints(grid) = length(first(RingGrids.get_londlatds(grid)))
_grid_npoints(spectral_grid::SpectralGrid) = spectral_grid.npoints
_grid_npoints(npoints::Integer) = npoints
