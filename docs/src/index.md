# SpeedyWeatherFTLE

SpeedyWeatherFTLE computes finite-time Lyapunov exponents (FTLEs) from
SpeedyWeather particle trajectories. The package can run a particle-tracking
simulation from prescribed velocity fields, compute positive- or negative-time
FTLE, reuse saved `ParticleTracker` NetCDF files, and convert FTLE arrays into
`RingGrids.Field` objects for plotting.

## What You Usually Need

Most workflows start with [`positive_FTLE`](@ref) or [`negative_FTLE`](@ref).
Pass zonal and meridional velocity fields on the same `RingGrids` grid and ask
for an [`FTLEResult`](@ref) when you want named fields plus metadata for
plotting and post-processing.

```@example quickstart
using CairoMakie
using Random
using RingGrids
using SpeedyWeatherFTLE

Random.seed!(42)

spatial_grid = FullGaussianGrid(8)
u = 25 * rand(spatial_grid)
v = 25 * rand(spatial_grid)

result = positive_FTLE(
    u,
    v;
    simulation_days = 1,
    dynamics = false,
    rint_hours = 6,
    return_result = true,
    time_indices = :nonzero,
)

size(result)
```

The result stores the selected FTLE integration horizons, the SpeedyWeather
spectral grid, the selected output times in hours, and run metadata:

```@example quickstart
result.direction, result.time_hours, result.dist_km
```

To plot the final selected output time, pass the result directly to
[`surface_plot`](@ref).

```@example quickstart
fig, ax, sp, cb = surface_plot(
    result;
    title = "Positive-time FTLE after $(result.time_hours[end]) hours",
    label = "FTLE [1/h]",
)

fig
```

For an interactive local integration-horizon plot, use [`slider_plot`](@ref).
In the static documentation build the slider is rendered but not interactive;
with GLMakie locally it is interactive.

```@example quickstart
fig, ax, sp, cb = slider_plot(
    result;
    title = "Positive-time FTLE",
    colorbar_label = "FTLE [1/h]",
)

fig
```

## Guide

- [Concepts and Data Layout](concepts.md): FTLE direction, units, particle
  layout, and `time_indices`.
- [Running Simulations](simulation.md): high-level simulation workflows with
  [`get_FTLE`](@ref), [`positive_FTLE`](@ref), and [`negative_FTLE`](@ref).
- [Particle Files](particle_files.md): saving and reusing SpeedyWeather
  `ParticleTracker` NetCDF output.
- [Plotting](plotting.md): converting arrays to fields and using
  [`surface_plot`](@ref), [`slider_plot`](@ref), [`animate_slider_plot`](@ref),
  and [`globe_plot`](@ref).
- [API Reference](api.md): generated reference documentation for exported
  functions and types.

## Development

From the repository root, instantiate the project once:

```julia
] instantiate
```

Run the package tests with:

```julia
] test
```

Build these docs locally with:

```bash
julia --project=docs docs/make.jl
```
