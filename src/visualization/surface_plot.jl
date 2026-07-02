"""
    surface_plot(field::RingGrids.Field; kwargs...)
    surface_plot(FTLE_grid::AbstractVector, grid_or_spectral_grid; kwargs...)
    surface_plot(FTLE_grid_time::AbstractMatrix, grid_or_spectral_grid; time_index = size(FTLE_grid_time, 2), kwargs...)
    surface_plot(FTLE_grid_time::AbstractMatrix, grid_or_spectral_grid; time_hours, time_hour = nothing, kwargs...)
    surface_plot(result::FTLEResult; time_index = size(result.ftle, 2), time_hour = nothing, kwargs...)

Plot one FTLE field on a geographic Makie axis.

`surface_plot` accepts an existing `RingGrids.Field`, a single-time FTLE
vector, the matrix returned by [`get_FTLE`](@ref), or an [`FTLEResult`](@ref).
For array inputs, pass either the `SpectralGrid` returned by the simulation API
or its spatial grid. FTLE array and result inputs label the colorbar as
`FTLE [1/h]` and use finite FTLE color limits by default; pass
`label = nothing` to suppress the label or `colorrange = nothing` for Makie
autoscaling.

# Keyword Arguments

- `lon = Vector(-180:180)`: interpolation longitudes for plotting.
- `lat = Vector(-90:90)`: interpolation latitudes for plotting.
- `shading = NoShading`: Makie surface shading option.
- `title = nothing`: optional plot title.
- `colormap = :viridis`: Makie colormap.
- `colorrange = nothing`: optional color limits for `Field` inputs; FTLE inputs default to finite-value extrema. Use [`ftle_colorrange`](@ref) for comparable plots.
- `colorbar = true`: add a colorbar.
- `label = nothing`: optional colorbar label for `Field` inputs; FTLE inputs default to `FTLE [1/h]`.
- `coastlines = true`: draw GeoMakie coastlines as a visual overlay.
- `coastline_color = :black`: coastline color.
- `coastline_linewidth = 1`: coastline line width.
- `figure_kwargs = (;)`: extra keyword arguments forwarded to `Figure`.
- `axis_kwargs = (;)`: extra keyword arguments forwarded to `GeoAxis`.
- `surface_kwargs = (;)`: extra keyword arguments forwarded to `surface!`.
- `colorbar_kwargs = (;)`: extra keyword arguments forwarded to `Colorbar`.
- `coastline_kwargs = (;)`: extra keyword arguments forwarded to `lines!`.
- `time_index = size(FTLE_grid_time, 2)`: selected FTLE column for matrix or
  result inputs.
- `time_hours = nothing`: saved integration horizons for matrix inputs when
  selecting with `time_hour`.
- `time_hour = nothing`: for matrix inputs with `time_hours`, or for
  `FTLEResult` inputs, select the saved integration horizon nearest this hour.

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
    field = ftle_field(FTLE_grid, grid_or_spectral_grid)
    plot_kwargs = (; kwargs...)
    if !(:label in keys(plot_kwargs))
        plot_kwargs = merge((; label=_FTLE_COLORBAR_LABEL), plot_kwargs)
    end
    if !(:colorrange in keys(plot_kwargs))
        plot_kwargs = merge((; colorrange=ftle_colorrange(FTLE_grid)), plot_kwargs)
    end
    return surface_plot(field; plot_kwargs...)
end

function surface_plot(
    FTLE_grid_time::AbstractMatrix,
    grid_or_spectral_grid;
    time_index::Union{Nothing,Integer}=nothing,
    time_hour::Union{Nothing,Real}=nothing,
    time_hours=nothing,
    kwargs...
    )
    resolved_index = _resolve_time_index(FTLE_grid_time; time_index, time_hour, time_hours)

    ftle_values = view(FTLE_grid_time, :, resolved_index)
    field = ftle_field(ftle_values, grid_or_spectral_grid)
    plot_kwargs = (; kwargs...)
    if !(:label in keys(plot_kwargs))
        plot_kwargs = merge((; label=_FTLE_COLORBAR_LABEL), plot_kwargs)
    end
    if !(:colorrange in keys(plot_kwargs))
        plot_kwargs = merge((; colorrange=ftle_colorrange(ftle_values)), plot_kwargs)
    end
    return surface_plot(field; plot_kwargs...)
end

function surface_plot(
    result::FTLEResult;
    time_index::Union{Nothing,Integer}=nothing,
    time_hour::Union{Nothing,Real}=nothing,
    kwargs...
    )
    spectral_grid = _require_spectral_grid(result, "surface_plot")
    return surface_plot(result.ftle, spectral_grid; time_index, time_hour, time_hours=result.time_hours, kwargs...)
end

export surface_plot
