"""
    SliderPlotHandle

Controls returned by [`slider_plot`](@ref) when `return_handle = true`.

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

function _plotting_extension_error(function_name)
    throw(ArgumentError(
        "$function_name requires a Makie backend and GeoMakie. Load them first, " *
        "for example `using CairoMakie, GeoMakie` or `using GLMakie, GeoMakie`.",
    ))
end

"""
    surface_plot(args...; kwargs...)

Plot one `RingGrids.Field`, FTLE array, or [`FTLEResult`](@ref) on a geographic
Makie axis.

Load GeoMakie and a Makie backend before calling this function, for example
`using CairoMakie, GeoMakie` for static output or `using GLMakie, GeoMakie` for
interactive windows. Plotting methods are loaded only after GeoMakie is
available, so compute-only users can install and load SpeedyWeatherFTLE without
the plotting stack.
"""
surface_plot(args...; kwargs...) = _plotting_extension_error("surface_plot")

"""
    slider_plot(args...; kwargs...)

Plot time-dependent fields or FTLE integration horizons with a Makie slider.

For FTLE arrays and [`FTLEResult`](@ref) inputs, the zero-duration sample is
skipped by default because FTLE is undefined at `t = 0`; signed nonzero
durations are supported. Load GeoMakie and a Makie backend before calling.
"""
slider_plot(args...; kwargs...) = _plotting_extension_error("slider_plot")

"""
    set_slider_time!(handle::SliderPlotHandle, time_hour)

Move a [`SliderPlotHandle`](@ref) to the saved time nearest to `time_hour`.
"""
set_slider_time!(args...; kwargs...) = _plotting_extension_error("set_slider_time!")

"""
    animate_slider_plot(path, args...; kwargs...)

Record an animation by advancing the same slider used by [`slider_plot`](@ref).
This requires GeoMakie and an active Makie backend.
"""
animate_slider_plot(args...; kwargs...) = _plotting_extension_error("animate_slider_plot")

"""
    globe_plot(args...; kwargs...)

Plot one `RingGrids.Field`, FTLE array, or [`FTLEResult`](@ref) on a GeoMakie
`GlobeAxis`. With GLMakie active, the returned globe can be rotated and zoomed;
with CairoMakie it renders as a static figure.
"""
globe_plot(args...; kwargs...) = _plotting_extension_error("globe_plot")

export SliderPlotHandle
export surface_plot
export slider_plot
export set_slider_time!
export animate_slider_plot
export globe_plot
