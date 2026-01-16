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
        lines!(ax, GeoMakie.coastlines(), color=:black)
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

export surface_plot
