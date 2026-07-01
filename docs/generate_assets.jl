using CairoMakie
using RingGrids
using SpeedyWeatherFTLE

assets_dir = joinpath(@__DIR__, "src", "assets")
mkpath(assets_dir)

spatial_grid = FullGaussianGrid(8)
londs, latds = RingGrids.get_londlatds(spatial_grid)

synthetic_ftle = @. 0.015 + 0.005 * sind(latds)^2
time_hours = [6.0, 12.0, 18.0]
FTLE_grid_time = hcat(synthetic_ftle, 1.5 .* synthetic_ftle, 2 .* synthetic_ftle)

animate_slider_plot(
    joinpath(assets_dir, "synthetic-ftle-slider.gif"),
    time_hours,
    FTLE_grid_time,
    spatial_grid;
    frames = eachindex(time_hours),
    framerate = 2,
    title = "Synthetic FTLE animation",
    colorbar_label = "FTLE [1/h]",
    coastlines = false,
)
