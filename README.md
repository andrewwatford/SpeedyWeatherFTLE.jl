# SpeedyWeatherFTLE

SpeedyWeatherFTLE computes finite-time Lyapunov exponent (FTLE) fields from
SpeedyWeather flow fields. Internally it uses particle advection, but the public
API is about the flow: call `FTLE`, get `RingGrids.Field` values, then analyze
or plot the fields.

The package stays narrow:

- compute forward- or backward-time FTLE with `FTLE(u, v; backwards = ...)`;
- return FTLE as `RingGrids.Field` time series plus selected integration times;
- compute finite-time stretching factors;
- plot FTLE fields with optional Makie/GeoMakie surface, slider, and globe
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

grid = FullGaussianGrid(4)
u = 0 .* rand(grid)
v = 0 .* rand(grid)

ftle, time_hours = FTLE(
    u,
    v;
    simulation_days = 0.25,
    rint_hours = 3,
    time_indices = :last,
)

final_field = ftle[:, end]
stretch = stretching_factor(ftle, time_hours)
```

Use `backwards = true` for backward-time FTLE:

```julia
backward_ftle, backward_time_hours = FTLE(
    u,
    v;
    backwards = true,
    simulation_days = 0.25,
    rint_hours = 3,
    time_indices = :last,
)
```

`FTLE` defaults to `radius = SpeedyWeatherFTLE.Re`, which is inherited from
SpeedyWeather. Pass a different radius in metres to change the spherical
distance scale used by the internal model and FTLE particle geometry.

Dynamic model workflows are out of scope for this package. `FTLE` initializes
the internal SpeedyWeather simulation with `dynamics = false`.

## Plotting

Load GeoMakie and a backend before using the plotting helpers:

```julia
using CairoMakie
using GeoMakie
using SpeedyWeatherFTLE

CairoMakie.activate!()

fig, ax, sp, cb = surface_plot(ftle; time_hours, coastlines = false)
handle = slider_plot(time_hours, ftle; return_handle = true, coastlines = false)
set_slider_time!(handle, last(time_hours))
fig_globe, ax_globe, sp_globe, cb_globe = globe_plot(ftle; time_hours, coastlines = false)
```

## Public API

- `FTLE`
- `stretching_factor`, `stretching_factor!`
- `surface_plot`, `slider_plot`, `globe_plot`, `animate_slider_plot`
- `SliderPlotHandle`, `set_slider_time!`, `shared_colorrange`
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
