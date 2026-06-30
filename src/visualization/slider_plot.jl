function _finite_absmax(data)
    max_val = 0.0
    found_finite = false
    for value in data
        if isfinite(value)
            max_val = max(max_val, abs(float(value)))
            found_finite = true
        end
    end
    return found_finite && max_val > 0 ? max_val : 1.0
end

"""
    slider_plot(times, field_ts::RingGrids.Field; kwargs...)
    slider_plot(times, FTLE_grid_time::AbstractMatrix, grid_or_spectral_grid; start_index = nothing, kwargs...)
    slider_plot(result::FTLEResult; kwargs...)

Plot an FTLE time series with a Makie slider.

`slider_plot` accepts a time-dependent `RingGrids.Field`, the matrix returned by
[`get_FTLE`](@ref), or an [`FTLEResult`](@ref). For FTLE arrays, the
zero-duration sample is skipped by default because FTLE is undefined at
`t = 0`; pass `start_index = 1` to include it.

# Keyword Arguments

- `lon = Vector(-180:180)`: interpolation longitudes for plotting.
- `lat = Vector(-90:90)`: interpolation latitudes for plotting.
- `shading = NoShading`: Makie surface shading option.
- `title = nothing`: optional plot title.
- `colormap = :viridis`: Makie colormap.
- `colorbar = true`: add a colorbar.
- `colorrange = nothing`: color limits. When omitted, finite values determine a symmetric range.
- `colorbar_label = nothing`: optional colorbar label.
- `coastlines = true`: draw GeoMakie coastlines.

# Returns

`fig, ax, sp, cb`, where `cb` is `nothing` when `colorbar = false`.
"""
function slider_plot(
    times::AbstractVector{<:Real},
    field_ts::Field;
    lon::Vector=Vector(-180:180),
    lat::Vector=Vector(-90:90),
    shading=NoShading, 
    title=nothing,
    colormap=:viridis,
    colorbar::Bool=true,
    colorrange=nothing,
    colorbar_label=nothing, 
    coastlines::Bool=true,
    )
    """
    Create a surface plot of a 2D field on a geographic axis. Currently supports only the default projection in GeoMakie.

    # Arguments
    - `field_ts::Field`: 2D field data to be plotted, in a time-series.
    - `duration::Float64`: Duration of the time series.
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
    # Number of time steps
    n_times = size(field_ts, 2)
    length(times) == n_times ||
        throw(DimensionMismatch("times has length $(length(times)), but field_ts has $n_times time steps"))
    # Set up longitude and latitude vectors
    lon_vec = vec([lo for lo in lon, la in lat])
    lat_vec = vec([la for lo in lon, la in lat])
    # Create figure and axis
    fig = Figure()
    if title !== nothing
        ax = GeoAxis(fig[1,1]; title=title)
    else
        ax = GeoAxis(fig[1,1])
    end

    # Construct the colorrange if not provided
    if colorrange === nothing
        max_val = _finite_absmax(field_ts)
        colorrange = (-max_val, max_val)
    end
    # And the colorbar if requested
    if colorbar
        if colorbar_label === nothing
            cb = Colorbar(fig[1, 2]; height=Relative(0.7), colorrange=colorrange)
        else
            cb = Colorbar(fig[1, 2]; label=colorbar_label, height=Relative(0.7), colorrange=colorrange)
        end
    else
        cb = nothing
    end

    # Set up the slider
    sg = SliderGrid(
        fig[2, 1],
        (label = "Time [h]", range = 1:n_times, format = i -> "$(times[Int(i)]) h", startvalue = 1))
    sl = sg.sliders[1]

    # Define the lifted field data for the surface plot
    field_data = lift(sl.value) do idx
        field = field_ts[:, idx]
        interpolate(lon_vec, lat_vec, field)
    end

    sp = surface!(ax, lon_vec, lat_vec, field_data; shading=shading, colormap=colormap, colorrange=colorrange)
    if coastlines
        lines!(ax, GeoMakie.coastlines(), color=:black, overdraw=true)
    end

    return fig, ax, sp, cb
end

function slider_plot(
    times::AbstractVector{<:Real},
    FTLE_grid_time::AbstractMatrix,
    grid_or_spectral_grid;
    start_index=nothing,
    kwargs...
    )
    """
    Create a slider plot directly from the `(FTLE_grid_time, spectral_grid,
    time_hours)` values returned by `get_FTLE`.

    By default the zero-duration sample is skipped, since its FTLE is undefined.
    Set `start_index=1` to include it.
    """
    length(times) == size(FTLE_grid_time, 2) ||
        throw(DimensionMismatch("times has length $(length(times)), but FTLE_grid_time has $(size(FTLE_grid_time, 2)) time steps"))

    if start_index === nothing
        start_index = something(findfirst(t -> isfinite(t) && t > 0, times), firstindex(times))
    end
    firstindex(times) <= start_index <= lastindex(times) ||
        throw(BoundsError(times, start_index))

    time_indices = start_index:lastindex(times)
    field_ts = ftle_field(view(FTLE_grid_time, :, time_indices), grid_or_spectral_grid)
    return slider_plot(times[time_indices], field_ts; kwargs...)
end

function slider_plot(
    result::FTLEResult;
    kwargs...
    )
    """
    Create a slider plot from an `FTLEResult`.
    """
    return slider_plot(result.time_hours, result.ftle, result.spectral_grid; kwargs...)
end

export slider_plot
