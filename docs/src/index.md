# SpeedyWeatherFTLE

SpeedyWeatherFTLE computes finite-time Lyapunov exponent (FTLE) fields from
SpeedyWeather flow fields. Users work with flow fields and returned
`RingGrids.Field` values; the particle advection used to estimate deformation is
an implementation detail.

## Install

SpeedyWeatherFTLE requires Julia 1.12 and the `mk/lyapunov2` SpeedyWeather
source branch because that branch provides the particle advection API used
internally.

```julia-repl
julia> using Pkg

julia> speedyweather_url = "https://github.com/SpeedyWeather/SpeedyWeather.jl";

julia> Pkg.add([
           PackageSpec(url = speedyweather_url, rev = "mk/lyapunov2", subdir = "LowerTriangularArrays"),
           PackageSpec(url = speedyweather_url, rev = "mk/lyapunov2", subdir = "RingGrids"),
           PackageSpec(url = speedyweather_url, rev = "mk/lyapunov2", subdir = "SpeedyTransforms"),
           PackageSpec(url = speedyweather_url, rev = "mk/lyapunov2", subdir = "SpeedyWeather"),
           PackageSpec(url = speedyweather_url, rev = "mk/lyapunov2", subdir = "SpeedyWeatherInternals"),
       ])

julia> Pkg.add(PackageSpec(url = "https://github.com/andrewwatford/SpeedyWeatherFTLE.jl"))
```

For plotting:

```julia-repl
julia> Pkg.add(["CairoMakie", "GeoMakie"])
```

## First Calculation

```@example first
using Logging
using RingGrids
using SpeedyWeatherFTLE

grid = FullGaussianGrid(4)
u = 0 .* rand(grid)
v = 0 .* rand(grid)

ftle, time_hours = with_logger(NullLogger()) do
    FTLE(
        u,
        v;
        simulation_days = 0.25,
        rint_hours = 3,
        particle_advection_every_n_time_steps = 1,
        time_indices = :last,
    )
end

(ftle isa Field, size(ftle), time_hours)
```

`ftle` is a `RingGrids.Field` with dimensions `(grid point, selected time)`.
The grid is carried by the field itself, so the matching time vector is the only
extra return value.

```@example first
final_field = ftle[:, end]
stretch = stretching_factor(ftle, time_hours)

(final_field isa Field, stretch isa Field)
```

Use `backwards = true` for backward-time FTLE.

## What To Read Next

- [Running FTLE](@ref) gives the full flow-field workflow and return values.
- [Plotting](@ref) shows `surface_plot`, `slider_plot`, and `globe_plot` with
  returned FTLE fields.
