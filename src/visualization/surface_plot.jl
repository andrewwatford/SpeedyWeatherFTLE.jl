"""
    surface_plot(field::RingGrids.Field; kwargs...)
    surface_plot(FTLE_grid::AbstractVector, grid_or_spectral_grid; kwargs...)
    surface_plot(FTLE_grid_time::AbstractMatrix, grid_or_spectral_grid; time_index = size(FTLE_grid_time, 2), kwargs...)
    surface_plot(result::FTLEResult; time_index = size(result.ftle, 2), kwargs...)

Plot one FTLE field on a geographic Makie axis.

`surface_plot` accepts an existing `RingGrids.Field`, a single-time FTLE
vector, the matrix returned by [`get_FTLE`](@ref), or an [`FTLEResult`](@ref).
For array inputs, pass either the `SpectralGrid` returned by the simulation API
or its spatial grid.

# Keyword Arguments

- `lon = Vector(-180:180)`: interpolation longitudes for plotting.
- `lat = Vector(-90:90)`: interpolation latitudes for plotting.
- `shading = NoShading`: Makie surface shading option.
- `title = nothing`: optional plot title.
- `colormap = :viridis`: Makie colormap.
- `colorbar = true`: add a colorbar.
- `label = nothing`: optional colorbar label.
- `coastlines = true`: draw GeoMakie coastlines.

# Returns

`fig, ax, sp, cb`, where `cb` is `nothing` when `colorbar = false`.
"""
function surface_plot(
    field::Field;
    lon::Vector=Vector(-180:180),
    lat::Vector=Vector(-90:90),
    shading=NoShading, 
    title=nothing,
    colormap=:viridis,
    colorbar::Bool=true,
    label=nothing, 
    coastlines::Bool=true,
    )
    """
    Create a surface plot of a 2D field on a geographic axis. Currently supports only the default projection in GeoMakie.

    # Arguments
    - `field::Field`: 2D field data to be plotted.
    - `lon::Vector`: Longitudes of the points to interpolate onto.
    - `lat::Vector`: Latitudes of the points to interpolate onto.
    - `shading`: Shading option for the surface plot (default: `NoShading`).
    - `title`: Title of the plot. (default: `nothing`).
    - `colormap`: Colormap to use for the surface plot. (default: `:viridis`).
    - `colorbar::Bool`: Whether to include a colorbar in the plot. (default: `true`).
    - `label`: Label for the colorbar. (default: `nothing`).
    - `coastlines::Bool`: Whether to add coastlines to the plot. (default: `true`).

    # Returns
    - `fig`: The figure object.
    - `ax`: The geographic axis object.
    - `sp`: The surface plot object.
    - `cb`: The colorbar object (if `colorbar` is `true`).
    """
    lon_vec = vec([lo for lo in lon, la in lat])
    lat_vec = vec([la for lo in lon, la in lat])
    field_data = interpolate(lon_vec, lat_vec, field)
    fig = Figure()
    if title !== nothing
        ax = GeoAxis(fig[1,1]; title=title)
    else
        ax = GeoAxis(fig[1,1])
    end
    sp = surface!(ax, lon_vec, lat_vec, field_data; shading=shading, colormap=colormap)
    if coastlines
        lines!(ax, GeoMakie.coastlines(), color=:black, overdraw=true)
    end
    if colorbar
        if label === nothing
            cb = Colorbar(fig[1, 2], sp; height=Relative(0.7))
        else
            cb = Colorbar(fig[1, 2], sp; label=label, height=Relative(0.7))
        end
    else
        cb = nothing
    end
    return fig, ax, sp, cb
end

function surface_plot(
    FTLE_grid::AbstractVector,
    grid_or_spectral_grid;
    kwargs...
    )
    """
    Create a surface plot directly from a single-time FTLE vector.

    `grid_or_spectral_grid` may be either the `SpectralGrid` returned by
    `get_FTLE` or its spatial grid.
    """
    field = ftle_field(FTLE_grid, grid_or_spectral_grid)
    return surface_plot(field; kwargs...)
end

function surface_plot(
    FTLE_grid_time::AbstractMatrix,
    grid_or_spectral_grid;
    time_index::Integer=size(FTLE_grid_time, 2),
    kwargs...
    )
    """
    Create a surface plot directly from the `FTLE_grid_time` matrix returned by
    `get_FTLE`.

    `grid_or_spectral_grid` may be either the `SpectralGrid` returned by
    `get_FTLE` or its spatial grid.
    """
    1 <= time_index <= size(FTLE_grid_time, 2) ||
        throw(BoundsError(FTLE_grid_time, (:, time_index)))

    field = ftle_field(view(FTLE_grid_time, :, time_index), grid_or_spectral_grid)
    return surface_plot(field; kwargs...)
end

function surface_plot(
    result::FTLEResult;
    time_index::Integer=size(result.ftle, 2),
    kwargs...
    )
    """
    Create a surface plot from an `FTLEResult`.
    """
    return surface_plot(result.ftle, result.spectral_grid; time_index, kwargs...)
end

export surface_plot
