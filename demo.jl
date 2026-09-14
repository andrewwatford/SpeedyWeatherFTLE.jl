# using Pkg

# speedyweather_url = "https://github.com/SpeedyWeather/SpeedyWeather.jl";

# Pkg.add([
#            PackageSpec(url = speedyweather_url, rev = "mk/lyapunov2", subdir = "LowerTriangularArrays"),
#            PackageSpec(url = speedyweather_url, rev = "mk/lyapunov2", subdir = "RingGrids"),
#            PackageSpec(url = speedyweather_url, rev = "mk/lyapunov2", subdir = "SpeedyTransforms"),
#            PackageSpec(url = speedyweather_url, rev = "mk/lyapunov2", subdir = "SpeedyWeather"),
#            PackageSpec(url = speedyweather_url, rev = "mk/lyapunov2", subdir = "SpeedyWeatherInternals"),
#        ])

# Pkg.add(PackageSpec(url = "https://github.com/andrewwatford/SpeedyWeatherFTLE.jl"))

# Pkg.add(["CairoMakie", "GeoMakie"])

using CairoMakie
using GeoMakie
using RingGrids
using SpeedyWeather
using SpeedyWeatherFTLE

spectral_grid = SpectralGrid(trunc=120, nlayers=1)
still_earth = Earth(spectral_grid, rotation=0)
initial_conditions = RandomVelocity(spectral_grid)
forcing = nothing
drag = nothing
model = BarotropicModel(spectral_grid; initial_conditions, planet=still_earth, forcing, drag)
simulation = initialize!(model)

plot_lon = -180:2:180
plot_lat = -90:2:90

for t in 1:100
    run!(simulation, period=Day(1))
    u = copy(simulation.variables.grid.u[:, 1])
    v = copy(simulation.variables.grid.v[:, 1])
    ftle, _ = FTLE(
        u,
        v;
        time_indices = :last,
    )
    fig, ax, sp, cb = surface_plot(
        ftle;
        coastlines = false,
        colorbar = false,
        colorrange = (0, 0.03)
    )
    hidedecorations!()
    save("figs/ftle_$(t).png", fig)
end

gif_cmd = `magick -delay 10 $(for i in $(seq 1 1 100); do echo figs/ftle_${i}.png; done) -loop 0 anim.gif`
run(gif_cmd)