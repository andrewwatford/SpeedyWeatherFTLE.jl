module SpeedyWeatherFTLEMakieExt

using GeoMakie
using RingGrids: Field, interpolate
using SpeedyWeatherFTLE

import SpeedyWeatherFTLE:
    SliderPlotHandle,
    animate_slider_plot,
    globe_plot,
    set_slider_time!,
    shared_colorrange,
    slider_plot,
    surface_plot

const SWFTLE = SpeedyWeatherFTLE
const Makie = GeoMakie.Makie

function _drop_keyword(kwargs::NamedTuple, key::Symbol)
    names = Tuple(name for name in keys(kwargs) if name != key)
    return NamedTuple{names}(map(name -> getproperty(kwargs, name), names))
end

function _normalize_slider_label_kwargs(kwargs)
    plot_kwargs = (; kwargs...)
    has_label = :label in keys(plot_kwargs)
    has_colorbar_label = :colorbar_label in keys(plot_kwargs)
    if has_label && has_colorbar_label
        throw(ArgumentError("pass either label or colorbar_label to slider_plot, not both"))
    elseif has_label
        colorbar_label = plot_kwargs.label
        plot_kwargs = _drop_keyword(plot_kwargs, :label)
        plot_kwargs = merge((; colorbar_label), plot_kwargs)
    end
    return plot_kwargs
end

function _resolve_colorrange(data, colorrange)
    colorrange === nothing && return shared_colorrange(data)
    colorrange === :auto && return shared_colorrange(data)
    colorrange === :symmetric && return shared_colorrange(data; symmetric=true)
    return colorrange
end

function _lonlat_vectors(lon, lat)
    lon_vec = vec([lo for lo in lon, la in lat])
    lat_vec = vec([la for lo in lon, la in lat])
    return lon_vec, lat_vec
end

function _field_on_lonlat(field::Field, lon, lat)
    lon_vec, lat_vec = _lonlat_vectors(lon, lat)
    return reshape(interpolate(lon_vec, lat_vec, field), length(lon), length(lat))
end

function _plot_tuple(handle::SliderPlotHandle)
    return handle.fig, handle.ax, handle.sp, handle.cb
end

function _selected_field(field::Field; time_index=nothing, time_hour=nothing, time_hours=nothing)
    if ndims(field) == 1
        time_hour === nothing || throw(ArgumentError("time_hour requires a time-dependent field"))
        (time_index === nothing || time_index == 1) || throw(BoundsError(field, (:, time_index)))
        return field
    end

    index = SWFTLE._resolve_time_index(field; time_index, time_hour, time_hours)
    return field[:, index]
end

"""
    surface_plot(field::RingGrids.Field; kwargs...)
    surface_plot(field_time::RingGrids.Field; time_index = nothing, time_hours = nothing, time_hour = nothing, kwargs...)

Plot one field or selected FTLE horizon on a geographic Makie axis.
"""
function surface_plot(
    field::Field;
    time_index::Union{Nothing,Integer}=nothing,
    time_hour::Union{Nothing,Real}=nothing,
    time_hours=nothing,
    lon::AbstractVector=Vector(-180:180),
    lat::AbstractVector=Vector(-90:90),
    shading=Makie.NoShading,
    title=nothing,
    colormap=:viridis,
    colorrange=:auto,
    colorbar::Bool=true,
    label=SWFTLE._FTLE_COLORBAR_LABEL,
    coastlines::Bool=true,
    coastline_color=:black,
    coastline_linewidth=1,
    figure_kwargs=NamedTuple(),
    axis_kwargs=NamedTuple(),
    surface_kwargs=NamedTuple(),
    colorbar_kwargs=NamedTuple(),
    coastline_kwargs=NamedTuple(),
)
    field = _selected_field(field; time_index, time_hour, time_hours)
    lon_vec, lat_vec = _lonlat_vectors(lon, lat)
    field_data = interpolate(lon_vec, lat_vec, field)

    fig = Makie.Figure(; figure_kwargs...)
    axis_attributes = title === nothing ? axis_kwargs : merge((; title), axis_kwargs)
    ax = GeoMakie.GeoAxis(fig[1, 1]; axis_attributes...)

    surface_attributes = merge((; shading, colormap), surface_kwargs)
    colorrange !== nothing && (surface_attributes = merge(surface_attributes, (; colorrange=_resolve_colorrange(field, colorrange))))
    sp = Makie.surface!(ax, lon_vec, lat_vec, field_data; surface_attributes...)

    if coastlines
        line_attributes = merge((; color=coastline_color, linewidth=coastline_linewidth, overdraw=true), coastline_kwargs)
        Makie.lines!(ax, GeoMakie.coastlines(); line_attributes...)
    end

    cb = if colorbar
        colorbar_attributes = merge((; height=Makie.Relative(0.7)), colorbar_kwargs)
        label === nothing ?
            Makie.Colorbar(fig[1, 2], sp; colorbar_attributes...) :
            Makie.Colorbar(fig[1, 2], sp; label, colorbar_attributes...)
    else
        nothing
    end

    return fig, ax, sp, cb
end

function _slider_plot_handle(
    times::AbstractVector{<:Real},
    field_ts::Field;
    lon::AbstractVector=Vector(-180:180),
    lat::AbstractVector=Vector(-90:90),
    shading=Makie.NoShading,
    title=nothing,
    colormap=:viridis,
    colorbar::Bool=true,
    colorrange=:auto,
    colorbar_label=nothing,
    coastlines::Bool=true,
    coastline_color=:black,
    coastline_linewidth=1,
    coastline_kwargs=NamedTuple(),
    slider_label="Time [h]",
    time_label::Bool=true,
    time_label_format=t -> "t = $(t) h",
    figure_kwargs=NamedTuple(),
    axis_kwargs=NamedTuple(),
    surface_kwargs=NamedTuple(),
    colorbar_kwargs=NamedTuple(),
)
    n_times = size(field_ts, 2)
    length(times) == n_times ||
        throw(DimensionMismatch("times has length $(length(times)), but field_ts has $n_times time steps"))

    lon_vec, lat_vec = _lonlat_vectors(lon, lat)
    fig = Makie.Figure(; figure_kwargs...)
    axis_attributes = title === nothing ? axis_kwargs : merge((; title), axis_kwargs)
    ax = GeoMakie.GeoAxis(fig[1, 1]; axis_attributes...)

    slider_layout = time_label ? Makie.GridLayout() : fig[2, 1]
    time_label && (fig[2, 1] = slider_layout)
    sg = Makie.SliderGrid(
        time_label ? slider_layout[2, 1] : slider_layout,
        (label=slider_label, range=1:n_times, format=i -> "$(times[Int(i)]) h", startvalue=1),
    )
    sl = sg.sliders[1]

    time_display = time_label ? Makie.Label(
        slider_layout[1, 1],
        Makie.lift(sl.value) do idx
            string(time_label_format(times[Int(idx)]))
        end;
        tellwidth=false,
    ) : nothing

    field_data = Makie.lift(sl.value) do idx
        interpolate(lon_vec, lat_vec, field_ts[:, idx])
    end

    surface_attributes = merge((; shading, colormap, colorrange=_resolve_colorrange(field_ts, colorrange)), surface_kwargs)
    sp = Makie.surface!(ax, lon_vec, lat_vec, field_data; surface_attributes...)

    if coastlines
        line_attributes = merge((; color=coastline_color, linewidth=coastline_linewidth, overdraw=true), coastline_kwargs)
        Makie.lines!(ax, GeoMakie.coastlines(); line_attributes...)
    end

    cb = if colorbar
        colorbar_attributes = merge((; height=Makie.Relative(0.7)), colorbar_kwargs)
        colorbar_label === nothing ?
            Makie.Colorbar(fig[1, 2], sp; colorbar_attributes...) :
            Makie.Colorbar(fig[1, 2], sp; label=colorbar_label, colorbar_attributes...)
    else
        nothing
    end

    return SliderPlotHandle(fig, ax, sp, cb, sg, sl, time_display, collect(times))
end

"""
    slider_plot(times, field_ts::RingGrids.Field; kwargs...)

Plot a time-dependent field or FTLE integration horizons with a Makie slider.
"""
function slider_plot(
    times::AbstractVector{<:Real},
    field_ts::Field;
    start_index=nothing,
    return_handle::Bool=false,
    kwargs...
)
    ndims(field_ts) == 2 ||
        throw(DimensionMismatch("slider_plot needs a time-dependent field with dimensions (grid point, time)"))
    length(times) == size(field_ts, 2) ||
        throw(DimensionMismatch("times has length $(length(times)), but field_ts has $(size(field_ts, 2)) columns"))

    if start_index === nothing
        start_index = findfirst(t -> isfinite(t) && !iszero(t), times)
        start_index === nothing &&
            throw(ArgumentError("no finite nonzero FTLE time is available; pass start_index=1 to include zero-duration data"))
    end
    firstindex(times) <= start_index <= lastindex(times) || throw(BoundsError(times, start_index))

    time_indices = start_index:lastindex(times)
    plot_kwargs = _normalize_slider_label_kwargs(kwargs)
    :colorbar_label in keys(plot_kwargs) || (plot_kwargs = merge((; colorbar_label=SWFTLE._FTLE_COLORBAR_LABEL), plot_kwargs))
    :slider_label in keys(plot_kwargs) || (plot_kwargs = merge((; slider_label=SWFTLE._FTLE_SLIDER_LABEL), plot_kwargs))
    :time_label_format in keys(plot_kwargs) || (plot_kwargs = merge((; time_label_format=t -> "Integration time = $(t) h"), plot_kwargs))

    handle = _slider_plot_handle(times[time_indices], field_ts[:, time_indices]; plot_kwargs...)
    return return_handle ? handle : _plot_tuple(handle)
end

function set_slider_time!(handle::SliderPlotHandle, time_hour::Real)
    Makie.set_close_to!(handle.slider, SWFTLE._nearest_time_index(handle.times, time_hour))
    return handle
end

function _record_slider_animation!(record_function, path, handle::SliderPlotHandle, frames; framerate, record_kwargs)
    record_function(handle.fig, path, frames; framerate, record_kwargs...) do frame_index
        Makie.set_close_to!(handle.slider, frame_index)
    end
    return path
end

function animate_slider_plot(
    path::AbstractString,
    times::AbstractVector{<:Real},
    field_ts::Field;
    start_index=nothing,
    framerate::Real=10,
    frames=nothing,
    record_kwargs=NamedTuple(),
    record_function=Makie.record,
    kwargs...
)
    handle = slider_plot(times, field_ts; start_index, return_handle=true, kwargs...)
    frame_indices = frames === nothing ? eachindex(handle.times) : frames
    return _record_slider_animation!(record_function, path, handle, frame_indices; framerate, record_kwargs)
end

"""
    globe_plot(field::RingGrids.Field; kwargs...)
    globe_plot(field_time::RingGrids.Field; time_index = nothing, time_hours = nothing, time_hour = nothing, kwargs...)

Plot one field or selected FTLE horizon on a GeoMakie `GlobeAxis`.
"""
function globe_plot(
    field::Field;
    time_index::Union{Nothing,Integer}=nothing,
    time_hour::Union{Nothing,Real}=nothing,
    time_hours=nothing,
    lon::AbstractVector=Vector(-180:180),
    lat::AbstractVector=Vector(-90:90),
    shading=Makie.NoShading,
    title=nothing,
    colormap=:viridis,
    colorrange=:auto,
    colorbar::Bool=true,
    label=SWFTLE._FTLE_COLORBAR_LABEL,
    coastlines::Bool=true,
    coastline_color=:black,
    coastline_linewidth=1,
    zlevel=10_000,
    coastline_zlevel=zlevel + 20_000,
    show_axis::Bool=false,
    camera_longlat=Makie.automatic,
    camera_altitude=Makie.automatic,
    figure_kwargs=NamedTuple(),
    axis_kwargs=NamedTuple(),
    surface_kwargs=NamedTuple(),
    colorbar_kwargs=NamedTuple(),
    coastline_kwargs=NamedTuple(),
)
    field = _selected_field(field; time_index, time_hour, time_hours)
    field_data = _field_on_lonlat(field, lon, lat)
    altitude = zeros(Float32, length(lon), length(lat))

    fig = Makie.Figure(; figure_kwargs...)
    axis_attributes = merge((; show_axis, title=title === nothing ? "" : title, camera_longlat, camera_altitude), axis_kwargs)
    ax = GeoMakie.GlobeAxis(fig[1, 1]; axis_attributes...)

    surface_attributes = merge((; color=field_data, shading, colormap, zlevel), surface_kwargs)
    colorrange !== nothing && (surface_attributes = merge(surface_attributes, (; colorrange=_resolve_colorrange(field, colorrange))))
    sp = Makie.surface!(ax, lon, lat, altitude; surface_attributes...)

    if coastlines
        line_attributes = merge(
            (; color=coastline_color, linewidth=coastline_linewidth, zlevel=coastline_zlevel, reset_limits=false),
            coastline_kwargs,
        )
        Makie.lines!(ax, GeoMakie.coastlines(); line_attributes...)
    end

    cb = if colorbar
        colorbar_attributes = merge((; height=Makie.Relative(0.7)), colorbar_kwargs)
        label === nothing ?
            Makie.Colorbar(fig[1, 2], sp; colorbar_attributes...) :
            Makie.Colorbar(fig[1, 2], sp; label, colorbar_attributes...)
    else
        nothing
    end

    return fig, ax, sp, cb
end

end
