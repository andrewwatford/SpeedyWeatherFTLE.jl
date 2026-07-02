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

Before running the examples in a fresh Julia project, add SpeedyWeatherFTLE and
the packages imported directly by the snippets:

```julia
]add https://github.com/andrewwatford/SpeedyWeatherFTLE.jl
]add RingGrids CairoMakie
```

For interactive GLMakie windows, add GLMakie separately:

```julia
]add GLMakie
```

Julia resolves `using PackageName` from the active project. If a script imports
`RingGrids`, `CairoMakie`, or `GLMakie` directly, add that package directly to
the active project even when it is also an indirect dependency. First
installation and precompilation can take several minutes with plotting
backends.

The repository development and docs environment uses `[sources]` entries to pin
SpeedyWeather monorepo packages to the `mk/lyapunov2` branch. A plain
`Pkg.add(url=...)` install resolves registered SpeedyWeather dependencies
instead. Use the local clone/develop setup below when you need to reproduce the
source-pinned development environment exactly.

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
[`surface_plot`](@ref). This quickstart uses random abstract velocity fields,
so coastlines are disabled; use them as visual context for geophysical fields.

```@example quickstart
fig, ax, sp, cb = surface_plot(
    result;
    title = "Positive-time FTLE after $(result.time_hours[end]) hours",
    label = "FTLE [1/h]",
    coastlines = false,
)

fig
```

For an interactive local integration-horizon plot, use [`slider_plot`](@ref).
In the static documentation build the slider is rendered but not interactive;
with GLMakie locally it is interactive. The slider sweeps the selected
integration durations from the same particle release; it is not an FTLE time
series.

```@example quickstart
fig, ax, sp, cb = slider_plot(
    result;
    title = "Positive-time FTLE",
    colorbar_label = "FTLE [1/h]",
    coastlines = false,
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

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
```

This develops the local checkout into the docs environment and instantiates the
docs dependencies.

For the package development environment itself:

```julia
] activate .
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

Local builds use non-pretty URLs so `docs/build/index.html` can be opened
directly from disk. Hosted CI builds use pretty URLs.
