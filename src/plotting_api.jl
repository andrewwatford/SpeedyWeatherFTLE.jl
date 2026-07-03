const _FTLE_COLORBAR_LABEL = "FTLE [1/h]"
const _FTLE_SLIDER_LABEL = "Integration time [h]"

function _plot_values(data)
    return data
end

"""
    shared_colorrange(data...; symmetric = false, pad = 0)

Return finite color limits spanning one or more arrays or `RingGrids.Field`
objects.
"""
function shared_colorrange(data...; symmetric::Bool=false, pad::Real=0)
    pad >= 0 || throw(ArgumentError("pad must be non-negative"))

    min_val = Inf
    max_val = -Inf
    found_finite = false
    for dataset in data
        for value in _plot_values(dataset)
            if isfinite(value)
                finite_value = float(value)
                min_val = min(min_val, finite_value)
                max_val = max(max_val, finite_value)
                found_finite = true
            end
        end
    end

    if !found_finite
        min_val, max_val = 0.0, 1.0
    elseif symmetric
        max_abs = max(abs(min_val), abs(max_val), eps(Float64))
        min_val, max_val = -max_abs, max_abs
    elseif min_val == max_val
        half_width = 0.05 * max(abs(min_val), 1.0)
        min_val, max_val = min_val - half_width, max_val + half_width
    end

    if pad > 0
        half_pad = pad * (max_val - min_val) / 2
        min_val -= half_pad
        max_val += half_pad
    end

    return (min_val, max_val)
end

function _plotting_extension_error(function_name)
    throw(ArgumentError(
        "$function_name requires GeoMakie and a Makie backend. Load them first, " *
        "for example `using CairoMakie, GeoMakie` or `using GLMakie, GeoMakie`.",
    ))
end

"""
    SliderPlotHandle

Controls returned by [`slider_plot`](@ref) when `return_handle = true`.
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

"""
    surface_plot(field::RingGrids.Field; kwargs...)
    surface_plot(field_ts::RingGrids.Field; time_index = nothing, time_hour = nothing, time_hours = nothing, kwargs...)

Plot one field or a selected column from a time-dependent `RingGrids.Field` on
a geographic Makie axis. For time-dependent fields, select a horizon with
`time_index` or nearest `time_hour`; by default the final column is shown.

Load GeoMakie plus a Makie backend before calling.
"""
surface_plot(args...; kwargs...) = _plotting_extension_error("surface_plot")

"""
    slider_plot(time_hours, field_ts::RingGrids.Field; start_index = nothing, return_handle = false, kwargs...)

Plot a time-dependent `RingGrids.Field` with a Makie slider. The field must have
dimensions `(grid point, time)` and `time_hours` must contain one value per
column. By default the slider starts at the first finite nonzero time.

Load GeoMakie plus a Makie backend before calling.
"""
slider_plot(args...; kwargs...) = _plotting_extension_error("slider_plot")

"""
    set_slider_time!(handle::SliderPlotHandle, time_hour)

Move a [`SliderPlotHandle`](@ref) to the saved time nearest `time_hour`.
"""
set_slider_time!(args...; kwargs...) = _plotting_extension_error("set_slider_time!")

"""
    animate_slider_plot(path, time_hours, field_ts::RingGrids.Field; start_index = nothing, framerate = 10, kwargs...)

Record an animation by advancing the same slider used by [`slider_plot`](@ref).
"""
animate_slider_plot(args...; kwargs...) = _plotting_extension_error("animate_slider_plot")

"""
    globe_plot(field::RingGrids.Field; kwargs...)
    globe_plot(field_ts::RingGrids.Field; time_index = nothing, time_hour = nothing, time_hours = nothing, kwargs...)

Plot one field or a selected column from a time-dependent `RingGrids.Field` on a
GeoMakie `GlobeAxis`. For time-dependent fields, select a horizon with
`time_index` or nearest `time_hour`; by default the final column is shown.
"""
globe_plot(args...; kwargs...) = _plotting_extension_error("globe_plot")
