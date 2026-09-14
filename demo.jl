using Pkg

speedyweather_url = "https://github.com/SpeedyWeather/SpeedyWeather.jl";

Pkg.add([
           PackageSpec(url = speedyweather_url, rev = "mk/lyapunov2", subdir = "LowerTriangularArrays"),
           PackageSpec(url = speedyweather_url, rev = "mk/lyapunov2", subdir = "RingGrids"),
           PackageSpec(url = speedyweather_url, rev = "mk/lyapunov2", subdir = "SpeedyTransforms"),
           PackageSpec(url = speedyweather_url, rev = "mk/lyapunov2", subdir = "SpeedyWeather"),
           PackageSpec(url = speedyweather_url, rev = "mk/lyapunov2", subdir = "SpeedyWeatherInternals"),
       ])

Pkg.add(PackageSpec(url = "https://github.com/andrewwatford/SpeedyWeatherFTLE.jl"))

Pkg.add(["GLMakie", "GeoMakie"])

using GLMakie
using GeoMakie
using Logging
using SpeedyWeather
using SpeedyWeatherFTLE

GLMakie.activate!()

simulation_period = Day(10)
ftle_period_days = 1

spectral_grid = SpectralGrid(trunc = 21, nlayers = 1)
orography = EarthOrography(spectral_grid)
initial_conditions = ZonalJet(spectral_grid)
model = with_logger(NullLogger()) do
    ShallowWaterModel(spectral_grid; orography, initial_conditions)
end
simulation = with_logger(NullLogger()) do
    initialize!(model)
end
with_logger(NullLogger()) do
    run!(simulation, steps = 0)
end

with_logger(NullLogger()) do
    run!(simulation, period = simulation_period)
end

final_u = copy(simulation.variables.grid.u[:, 1])
final_v = copy(simulation.variables.grid.v[:, 1])

plot_lon = -180:2:180
plot_lat = -90:2:90

final_ftle, final_time_hours = with_logger(NullLogger()) do
    FTLE(
        final_u,
        final_v;
        simulation_days = ftle_period_days,
        rint_hours = 6,
        particle_advection_every_n_time_steps = 4,
        time_indices = :last,
    )
end

fig, ax, sp, cb = globe_plot(
    final_ftle;
    coastlines = false,
    colorbar = true,
)