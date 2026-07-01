"""
    ftle_field(FTLE_grid_time, grid_or_spectral_grid)
    ftle_field(FTLE_grid, grid_or_spectral_grid)
    ftle_field(result::FTLEResult; time_indices = :)

Convert FTLE arrays into `RingGrids.Field` objects.

Use the matrix form for a full `(grid point, integration horizon)` FTLE output,
the vector form for a single output time, or pass an [`FTLEResult`](@ref)
directly. The grid argument may be either the `SpectralGrid` returned by
[`get_FTLE`](@ref) or its spatial grid.
"""
function ftle_field(FTLE_grid_time::AbstractMatrix, grid_or_spectral_grid)
    """
    Convert the matrix returned by `get_FTLE` into a RingGrids `Field`.

    `grid_or_spectral_grid` may be either the `SpectralGrid` returned by
    `get_FTLE` or its spatial grid.
    """
    npoints = _grid_npoints(grid_or_spectral_grid)
    size(FTLE_grid_time, 1) == npoints ||
        throw(DimensionMismatch("FTLE_grid_time has $(size(FTLE_grid_time, 1)) rows, but the grid has $npoints points"))

    return Field(FTLE_grid_time, _spatial_grid(grid_or_spectral_grid))
end

function ftle_field(FTLE_grid::AbstractVector, grid_or_spectral_grid)
    """
    Convert a single-time FTLE vector into a RingGrids `Field`.

    `grid_or_spectral_grid` may be either the `SpectralGrid` returned by
    `get_FTLE` or its spatial grid.
    """
    npoints = _grid_npoints(grid_or_spectral_grid)
    length(FTLE_grid) == npoints ||
        throw(DimensionMismatch("FTLE_grid has length $(length(FTLE_grid)), but the grid has $npoints points"))

    return Field(FTLE_grid, _spatial_grid(grid_or_spectral_grid))
end

export ftle_field
