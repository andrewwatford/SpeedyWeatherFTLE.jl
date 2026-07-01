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

function _finite_extrema(data)
    min_val = Inf
    max_val = -Inf
    found_finite = false
    for value in data
        if isfinite(value)
            finite_value = float(value)
            min_val = min(min_val, finite_value)
            max_val = max(max_val, finite_value)
            found_finite = true
        end
    end

    if !found_finite
        return (0.0, 1.0)
    elseif min_val == max_val
        half_width = 0.05 * max(abs(min_val), 1.0)
        return (min_val - half_width, max_val + half_width)
    else
        return (min_val, max_val)
    end
end

function _resolve_colorrange(data, colorrange)
    if colorrange === nothing || colorrange === :auto
        return _finite_extrema(data)
    elseif colorrange === :symmetric
        max_val = _finite_absmax(data)
        return (-max_val, max_val)
    else
        return colorrange
    end
end

"""
    SliderPlotHandle

Internal controls returned by [`slider_plot`](@ref) when
`return_handle = true`.

The first four fields match the normal `slider_plot` return values:
`fig, ax, sp, cb`. The remaining fields expose the `SliderGrid`, the active
slider, the optional live time label, and the plotted time values so helper
functions such as [`set_slider_time!`](@ref) and [`animate_slider_plot`](@ref)
can drive the same slider plot.
"""
struct SliderPlotHandle{F, A, S, C, G, L, D, T}
    fig::F
    ax::A
    sp::S
    cb::C
    slidergrid::G
    slider::L
    time_label::D
    times::T
end

_plot_tuple(handle::SliderPlotHandle) = (handle.fig, handle.ax, handle.sp, handle.cb)

function _slider_plot_handle(
    times::AbstractVector{<:Real},
    field_ts::Field;
    lon::Vector=Vector(-180:180),
    lat::Vector=Vector(-90:90),
    shading=NoShading, 
    title=nothing,
    colormap=:viridis,
    colorbar::Bool=true,
    colorrange=:auto,
    colorbar_label=nothing, 
    coastlines::Bool=true,
    time_label::Bool=true,
    time_label_format=t -> "t = $(t) h",
    )
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

    # Construct the colorrange if not provided.
    resolved_colorrange = _resolve_colorrange(field_ts, colorrange)

    # Set up the slider
    slider_layout = time_label ? GridLayout() : fig[2, 1]
    if time_label
        fig[2, 1] = slider_layout
    end
    sg = SliderGrid(
        time_label ? slider_layout[2, 1] : slider_layout,
        (label = "Time [h]", range = 1:n_times, format = i -> "$(times[Int(i)]) h", startvalue = 1))
    sl = sg.sliders[1]

    time_display = if time_label
        Label(
            slider_layout[1, 1],
            lift(sl.value) do idx
                string(time_label_format(times[Int(idx)]))
            end;
            tellwidth=false,
        )
    else
        nothing
    end

    # Define the lifted field data for the surface plot
    field_data = lift(sl.value) do idx
        field = field_ts[:, idx]
        interpolate(lon_vec, lat_vec, field)
    end

    sp = surface!(ax, lon_vec, lat_vec, field_data; shading=shading, colormap=colormap, colorrange=resolved_colorrange)
    if coastlines
        lines!(ax, GeoMakie.coastlines(), color=:black, overdraw=true)
    end

    if colorbar
        if colorbar_label === nothing
            cb = Colorbar(fig[1, 2], sp; height=Relative(0.7))
        else
            cb = Colorbar(fig[1, 2], sp; label=colorbar_label, height=Relative(0.7))
        end
    else
        cb = nothing
    end

    return SliderPlotHandle(fig, ax, sp, cb, sg, sl, time_display, collect(times))
end

"""
    slider_plot(times, field_ts::RingGrids.Field; kwargs...)
    slider_plot(times, FTLE_grid_time::AbstractMatrix, grid_or_spectral_grid; start_index = nothing, kwargs...)
    slider_plot(result::FTLEResult; kwargs...)

Plot FTLE values for different integration horizons with a Makie slider.

`slider_plot` accepts a time-dependent `RingGrids.Field`, the matrix returned by
[`get_FTLE`](@ref), or an [`FTLEResult`](@ref). For FTLE arrays, the
zero-duration sample is skipped by default because FTLE is undefined at
`t = 0`; pass `start_index = 1` to include it. For FTLE outputs, the slider
shows different integration durations from the same particle release, not a
time series of independent instantaneous FTLE fields.

# Keyword Arguments

- `lon = Vector(-180:180)`: interpolation longitudes for plotting.
- `lat = Vector(-90:90)`: interpolation latitudes for plotting.
- `shading = NoShading`: Makie surface shading option.
- `title = nothing`: optional plot title.
- `colormap = :viridis`: Makie colormap.
- `colorbar = true`: add a colorbar.
- `colorrange = :auto`: color limits. Use `:auto` or `nothing` for finite-value extrema, `:symmetric` for symmetric limits, or pass explicit limits.
- `colorbar_label = nothing`: optional colorbar label.
- `coastlines = true`: draw GeoMakie coastlines.
- `time_label = true`: show a live label above the slider with the active time.
- `time_label_format = t -> "t = \$(t) h"`: format the live time label.
- `return_handle = false`: return a [`SliderPlotHandle`](@ref) with the slider
  controls instead of the usual four-value tuple.

# Returns

By default, returns `fig, ax, sp, cb`, where `cb` is `nothing` when
`colorbar = false`. With `return_handle = true`, returns a
[`SliderPlotHandle`](@ref).
"""
function slider_plot(
    times::AbstractVector{<:Real},
    field_ts::Field;
    return_handle::Bool=false,
    kwargs...
    )
    handle = _slider_plot_handle(times, field_ts; kwargs...)
    return return_handle ? handle : _plot_tuple(handle)
end

function slider_plot(
    times::AbstractVector{<:Real},
    FTLE_grid_time::AbstractMatrix,
    grid_or_spectral_grid;
    start_index=nothing,
    return_handle::Bool=false,
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
    return slider_plot(times[time_indices], field_ts; return_handle, kwargs...)
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

function _nearest_time_index(times, time_hour)
    isempty(times) && throw(ArgumentError("slider handle has no times"))
    nearest_index = firstindex(times)
    nearest_distance = abs(times[nearest_index] - time_hour)
    for index in Iterators.drop(eachindex(times), 1)
        distance = abs(times[index] - time_hour)
        if distance < nearest_distance
            nearest_index = index
            nearest_distance = distance
        end
    end
    return nearest_index
end

"""
    set_slider_time!(handle::SliderPlotHandle, time_hour)

Move a slider plot to the saved time nearest to `time_hour`.

Returns `handle`, so calls can be chained with display or animation code.
"""
function set_slider_time!(handle::SliderPlotHandle, time_hour::Real)
    time_index = _nearest_time_index(handle.times, time_hour)
    Makie.set_close_to!(handle.slider, time_index)
    return handle
end

function _record_slider_animation!(record_function, path, handle::SliderPlotHandle, frames; framerate, record_kwargs)
    record_function(handle.fig, path, frames; framerate, record_kwargs...) do frame_index
        Makie.set_close_to!(handle.slider, frame_index)
    end
    return path
end

"""
    animate_slider_plot(path, times, field_ts::RingGrids.Field; kwargs...)
    animate_slider_plot(path, times, FTLE_grid_time, grid_or_spectral_grid; start_index = nothing, kwargs...)
    animate_slider_plot(path, result::FTLEResult; kwargs...)

Record an animation by advancing the same slider used by [`slider_plot`](@ref).

`path` should end in a Makie-supported animation extension such as `.mp4` or
`.gif`. All plotting keyword arguments accepted by [`slider_plot`](@ref) are
forwarded. Use `framerate` to control the output speed. `record_kwargs` are
forwarded to `Makie.record`.

This function works with the active Makie backend. In CI and documentation
builds, `CairoMakie` can record a non-interactive animation. For local
interactive exploration, activate `GLMakie` before calling `slider_plot` or
`animate_slider_plot`.

# Returns

The animation `path`.
"""
function animate_slider_plot(
    path::AbstractString,
    times::AbstractVector{<:Real},
    field_ts::Field;
    framerate::Real=10,
    frames=nothing,
    record_kwargs=NamedTuple(),
    record_function=Makie.record,
    kwargs...
    )
    handle = slider_plot(times, field_ts; return_handle=true, kwargs...)
    frame_indices = frames === nothing ? eachindex(handle.times) : frames
    return _record_slider_animation!(record_function, path, handle, frame_indices; framerate, record_kwargs)
end

function animate_slider_plot(
    path::AbstractString,
    times::AbstractVector{<:Real},
    FTLE_grid_time::AbstractMatrix,
    grid_or_spectral_grid;
    start_index=nothing,
    framerate::Real=10,
    frames=nothing,
    record_kwargs=NamedTuple(),
    record_function=Makie.record,
    kwargs...
    )
    handle = slider_plot(
        times,
        FTLE_grid_time,
        grid_or_spectral_grid;
        start_index,
        return_handle=true,
        kwargs...,
    )
    frame_indices = frames === nothing ? eachindex(handle.times) : frames
    return _record_slider_animation!(record_function, path, handle, frame_indices; framerate, record_kwargs)
end

function animate_slider_plot(
    path::AbstractString,
    result::FTLEResult;
    kwargs...
    )
    return animate_slider_plot(path, result.time_hours, result.ftle, result.spectral_grid; kwargs...)
end

export slider_plot
export SliderPlotHandle
export set_slider_time!
export animate_slider_plot
