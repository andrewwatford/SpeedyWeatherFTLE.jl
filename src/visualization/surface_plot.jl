"""
    surface_plot(field::RingGrids.Field; kwargs...)
    surface_plot(FTLE_grid::AbstractVector, grid_or_spectral_grid; kwargs...)
    surface_plot(FTLE_grid_time::AbstractMatrix, grid_or_spectral_grid; time_index = size(FTLE_grid_time, 2), kwargs...)
    surface_plot(result::FTLEResult; time_index = size(result.ftle, 2), kwargs...)

Plot one FTLE field on a geographic Makie axis.

`surface_plot` accepts an existing `RingGrids.Field`, a single-time FTLE
vector, the matrix returned by [`get_FTLE`](@ref), or an [`FTLEResult`](@ref).
For array inputs, pass either the `SpectralGrid` returned by the simulation API
or its spatial grid. FTLE array and result inputs label the colorbar as
`FTLE [1/h]` by default; pass `label = nothing` to suppress it.

# Keyword Arguments

- `lon = Vector(-180:180)`: interpolation longitudes for plotting.
- `lat = Vector(-90:90)`: interpolation latitudes for plotting.
- `shading = NoShading`: Makie surface shading option.
- `title = nothing`: optional plot title.
- `colormap = :viridis`: Makie colormap.
- `colorrange = nothing`: optional color limits. Use [`ftle_colorrange`](@ref) for comparable plots.
- `colorbar = true`: add a colorbar.
- `label = nothing`: optional colorbar label for `Field` inputs; FTLE inputs default to `FTLE [1/h]`.
- `coastlines = true`: draw GeoMakie coastlines.
- `coastline_color = :black`: coastline color.
- `coastline_linewidth = 1`: coastline line width.
- `figure_kwargs = (;)`: extra keyword arguments forwarded to `Figure`.
- `axis_kwargs = (;)`: extra keyword arguments forwarded to `GeoAxis`.
- `surface_kwargs = (;)`: extra keyword arguments forwarded to `surface!`.
- `colorbar_kwargs = (;)`: extra keyword arguments forwarded to `Colorbar`.
- `coastline_kwargs = (;)`: extra keyword arguments forwarded to `lines!`.

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
    colorrange=nothing,
    colorbar::Bool=true,
    label=nothing, 
    coastlines::Bool=true,
    coastline_color=:black,
    coastline_linewidth=1,
    figure_kwargs=NamedTuple(),
    axis_kwargs=NamedTuple(),
    surface_kwargs=NamedTuple(),
    colorbar_kwargs=NamedTuple(),
    coastline_kwargs=NamedTuple(),
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
    - `colorrange`: Optional color limits.
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
    fig = Figure(; figure_kwargs...)
    axis_attributes = title === nothing ? axis_kwargs : merge((; title), axis_kwargs)
    ax = GeoAxis(fig[1,1]; axis_attributes...)

    surface_attributes = merge((; shading, colormap), surface_kwargs)
    if colorrange !== nothing
        surface_attributes = merge(surface_attributes, (; colorrange))
    end
    sp = surface!(ax, lon_vec, lat_vec, field_data; surface_attributes...)
    if coastlines
        line_attributes = merge(
            (; color=coastline_color, linewidth=coastline_linewidth, overdraw=true),
            coastline_kwargs,
        )
        lines!(ax, GeoMakie.coastlines(); line_attributes...)
    end
    if colorbar
        colorbar_attributes = merge((; height=Relative(0.7)), colorbar_kwargs)
        if label === nothing
            cb = Colorbar(fig[1, 2], sp; colorbar_attributes...)
        else
            cb = Colorbar(fig[1, 2], sp; label, colorbar_attributes...)
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
    if :label in keys(kwargs)
        return surface_plot(field; kwargs...)
    else
        return surface_plot(field; label=_FTLE_COLORBAR_LABEL, kwargs...)
    end
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
    if :label in keys(kwargs)
        return surface_plot(field; kwargs...)
    else
        return surface_plot(field; label=_FTLE_COLORBAR_LABEL, kwargs...)
    end
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
