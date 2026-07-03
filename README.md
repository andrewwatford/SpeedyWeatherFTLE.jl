# SpeedyWeatherFTLE

SpeedyWeatherFTLE computes finite-time Lyapunov exponent (FTLE) fields from
SpeedyWeather flow fields. Internally it uses particle advection, but the public
API is about the flow: call `FTLE`, get FTLE arrays or an `FTLEResult`, then
analyze or plot the result.

The package stays narrow:

- compute forward- or backward-time FTLE with `FTLE(u, v; backwards = ...)`;
- keep grid, time, and run metadata in `FTLEResult`;
- convert FTLE arrays to `RingGrids.Field`;
- compute finite-time stretching factors;
- plot FTLE results with optional Makie/GeoMakie surface, slider, and globe
  helpers.

The plotting layer is optional. Compute-only users do not load Makie.

## Installation

Requires Julia 1.12 or newer.

SpeedyWeatherFTLE currently needs particle-advection APIs from the
`mk/lyapunov2` branch of the SpeedyWeather monorepo:

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

For plotting, add GeoMakie plus a Makie backend:

```julia
Pkg.add(["CairoMakie", "GeoMakie"])
```

## Basic Usage

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

`FTLE` returns `(ftle, spectral_grid, time_hours)` by default. With
`return_result = true`, it returns an `FTLEResult`.

```julia
final_values = final_ftle(result)
final_field = final_ftle_field(result)
field_at_12h = ftle_field(result; time_hour = 12)
stretch = stretching_factor(result)
```

The supplied-field wrapper supports frozen prescribed flow fields. Leave
`dynamics = false`, or build a dedicated SpeedyWeather simulation for dynamic
model workflows.

## Plotting

Load GeoMakie and a backend before using the plotting helpers:

```julia
using CairoMakie
using GeoMakie
using SpeedyWeatherFTLE

CairoMakie.activate!()

fig, ax, sp, cb = surface_plot(result; coastlines = false)
handle = slider_plot(result; return_handle = true, coastlines = false)
set_slider_time!(handle, 12)
fig_globe, ax_globe, sp_globe, cb_globe = globe_plot(result; coastlines = false)
```

`surface_plot`, `slider_plot`, and `globe_plot` accept `FTLEResult` instances
directly. `surface_plot` and `globe_plot` also accept FTLE vectors or matrices
when passed with the corresponding grid.

## Public API

- `FTLE`
- `FTLEResult`, `final_ftle`, `final_ftle_field`, `ftle_field`
- `stretching_factor`, `stretching_factor!`
- `surface_plot`, `slider_plot`, `globe_plot`, `animate_slider_plot`
- `SliderPlotHandle`, `set_slider_time!`, `ftle_colorrange`
- `Re`

## Development

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
Pkg.test()
```

Build the documentation locally with:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```
