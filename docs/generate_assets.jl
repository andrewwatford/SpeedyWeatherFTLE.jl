using CairoMakie
using RingGrids
using SpeedyWeatherFTLE

assets_dir = joinpath(@__DIR__, "src", "assets")
mkpath(assets_dir)

spatial_grid = FullGaussianGrid(16)
londs, latds = RingGrids.get_londlatds(spatial_grid)

function turbulent_ftle(londs, latds, hour)
    phase = hour / 6
    meander = @. 32 + 8 * sind(2londs + 18phase)
    jet = @. exp(-((latds - meander) / 11)^2)
    shear_ridges = @. 0.010 * jet * (1 + 0.35 * cosd(5londs - 24phase))

    eddies = zeros(Float64, length(londs))
    for (lon0, lat0, width_lon, width_lat, spin) in (
        (-150, -24, 20, 12, 1),
        (-88, 18, 16, 10, -1),
        (-22, -8, 18, 14, 1),
        (44, 34, 22, 12, -1),
        (118, -30, 20, 11, 1),
        (158, 10, 15, 9, -1),
    )
        dlon = @. mod(londs - lon0 + 180, 360) - 180
        dlat = @. latds - lat0
        r2 = @. (dlon / width_lon)^2 + (dlat / width_lat)^2
        eddies .+= @. 0.007 * exp(-r2) * (1 + 0.45 * cosd(4dlon + 3dlat - 30spin * phase))
    end

    filaments = @. 0.0045 * abs(sind(3londs + 2latds + 22phase)) * exp(-abs(latds) / 65)
    growth = 1 - exp(-hour / 24)
    return @. 0.002 + growth * (abs(shear_ridges) + abs(eddies) + filaments)
end

time_hours = collect(6.0:6.0:72.0)
FTLE_grid_time = hcat((turbulent_ftle(londs, latds, hour) for hour in time_hours)...)

animate_slider_plot(
    joinpath(assets_dir, "synthetic-ftle-slider.gif"),
    time_hours,
    FTLE_grid_time,
    spatial_grid;
    frames = eachindex(time_hours),
    framerate = 6,
    title = "Synthetic turbulent FTLE horizons",
    colormap = :magma,
    coastline_color = :white,
    coastline_linewidth = 1.2,
)
