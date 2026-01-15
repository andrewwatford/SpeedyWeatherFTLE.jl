using SpeedyWeather, RingGrids
using GeoMakie, GLMakie

function slider_plot(
    times::Vector{Float64},
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
        max_val = maximum(abs.(field_ts))
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
        (label = "Time [h]", range = 1:n_times, format = i -> "$(times[i]) h", startvalue = 1))
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

export slider_plot