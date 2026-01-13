# 2D Turbulence

```@example
using SpeedyWeather, NCDatasets, GeoMakie, CairoMakie

# Set up the spectral grid
spectral_grid = SpectralGrid(trunc=63, nlayers=1)
# Initialize the barotropic model with random initial conditions
initial_conditions = RandomVelocity(spectral_grid)
# Initialize the model
model = BarotropicModel(spectral_grid; initial_conditions)
simulation = initialize!(model)
# Run the simulation for 20 days
run!(simulation, period=Day(20), output=true)

run_folder = model.output.run_folder
filename = model.output.filename
ds = NCDataset(joinpath(run_folder, filename))
lons = ds["lon"][:]
lats = ds["lat"][:]
vor = ds["vor"][:, :, 1, end]

fig = Figure()
ax = GeoAxis(fig[1, 1], title = "2D Turbulence: Vorticity after 20 days")
sp = surface!(ax, lons, lats, vor; shading = NoShading)
cb = Colorbar(fig[1, 2], sp; label = "ω [1/s]", width = 18, height = Relative(0.7))

fig
```
