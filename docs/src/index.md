# SpeedyWeatherFTLE

SpeedyWeatherFTLE computes finite-time Lyapunov exponent (FTLE) fields from
SpeedyWeather flow fields. Users work with flow fields and FTLE results; the
particle advection used to estimate the deformation is an implementation
detail.

## Install

SpeedyWeatherFTLE requires Julia 1.12 and the `mk/lyapunov2` SpeedyWeather
source branch because that branch provides the particle advection API used
internally.

```julia
using Pkg

speedyweather_url = "https://github.com/SpeedyWeather/SpeedyWeather.jl"
Pkg.add([
    PackageSpec(url = speedyweather_url, rev = "mk/lyapunov2", subdir = "LowerTriangularArrays"),
    PackageSpec(url = speedyweather_url, rev = "mk/lyapunov2", subdir = "RingGrids"),
    PackageSpec(url = speedyweather_url, rev = "mk/lyapunov2", subdir = "SpeedyTransforms"),
    PackageSpec(url = speedyweather_url, rev = "mk/lyapunov2", subdir = "SpeedyWeather"),
    PackageSpec(url = speedyweather_url, rev = "mk/lyapunov2", subdir = "SpeedyWeatherInternals"),
])
Pkg.add(PackageSpec(url = "https://github.com/andrewwatford/SpeedyWeatherFTLE.jl"))
```

For plotting:

```julia
Pkg.add(["CairoMakie", "GeoMakie"])
```

Use `CairoMakie` for static figures and documentation-style output. Use
`GLMakie` instead of `CairoMakie` for interactive local windows.

## First Calculation

```julia
using RingGrids
using SpeedyWeatherFTLE

grid = FullGaussianGrid(8)
u = 25 .* rand(grid)
v = 25 .* rand(grid)

result = FTLE(
    u,
    v;
    simulation_days = 1,
    rint_hours = 6,
    return_result = true,
    time_indices = :nonzero,
)
```

`result` is an [`FTLEResult`](@ref). It stores the FTLE matrix, the
`SpectralGrid`, selected integration horizons in hours, and run metadata.

```julia
final_values = final_ftle(result)
final_field = final_ftle_field(result)
field_at_12h = ftle_field(result; time_hour = 12)
stretch = stretching_factor(result)
```

Use `backwards = true` for backward-time FTLE:

```julia
backward = FTLE(
    u,
    v;
    backwards = true,
    simulation_days = 1,
    rint_hours = 6,
    return_result = true,
    time_indices = :last,
)
```

## What To Read Next

- [Running FTLE](@ref) gives the full flow-field workflow and return values.
- [Plotting](@ref) shows `surface_plot`, `slider_plot`, and `globe_plot` with
  direct [`FTLEResult`](@ref) dispatch.
