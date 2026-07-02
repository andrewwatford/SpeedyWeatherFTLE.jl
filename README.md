# SpeedyWeatherFTLE

[![Build Status](https://github.com/andrewwatford/SpeedyWeatherFTLE.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/andrewwatford/SpeedyWeatherFTLE.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![docs](https://img.shields.io/badge/documentation-latest_release-blue.svg)](https://andrewwatford.github.io/SpeedyWeatherFTLE.jl/)

SpeedyWeatherFTLE computes finite-time Lyapunov exponents (FTLEs) from
SpeedyWeather particle trajectories. It can run a SpeedyWeather particle
tracking simulation from prescribed velocity fields, compute positive- or
negative-time FTLE, reuse saved `ParticleTracker` NetCDF files, and convert
FTLE arrays to RingGrids fields for plotting.

## Installation

SpeedyWeatherFTLE currently depends on particle-tracking features from the
`mk/lyapunov2` branch of the SpeedyWeather monorepo. Install those source
dependencies first, then install SpeedyWeatherFTLE:

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
Pkg.add(["CairoMakie"])
```

Then load the packages in Julia with:

```julia
using CairoMakie
using RingGrids
using SpeedyWeatherFTLE
```

Julia resolves imports from the active project. Even when a package is an
indirect dependency of SpeedyWeatherFTLE, any package you `using` directly in a
script or notebook should be a direct dependency of that active project. A
plain `Pkg.add(url="https://github.com/andrewwatford/SpeedyWeatherFTLE.jl")`
without the `mk/lyapunov2` SpeedyWeather source dependencies is not a supported
install path.

For local interactive Makie windows and rotatable globes, add GLMakie
separately:

```julia
]add GLMakie
```

The first install and precompile can take several minutes, especially when
plotting backends are included.

## Documentation

The Documenter site in `docs/src` is the best onboarding path. It includes:

- concepts and array layout for FTLE post-processing;
- complete examples for positive- and negative-time FTLE;
- saved particle-file workflows;
- plotting examples with `surface_plot`, `slider_plot`, `animate_slider_plot`,
  and `globe_plot`;
- an API reference generated from the package docstrings.

Build it locally from the repository root with:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

Then open `docs/build/index.html` in your browser. If you prefer to browse with
a local web server, run:

```bash
python3 -m http.server --directory docs/build 8000
```

The checked-out docs and development environments use this repository's
`[sources]` entries for the SpeedyWeather monorepo packages on the
`mk/lyapunov2` branch. Fresh user installs must add those source dependencies
explicitly, as shown above, because Julia does not apply `[sources]` entries
transitively when this package is added as a dependency.

## Examples

The richest examples live in the documentation:

- `docs/src/simulation.md` shows SpeedyWeather particle-tracking workflows,
  including a summer-like versus winter-like meandering jet comparison and an
  evolving-flow versus frozen-flow recipe.
- `docs/src/plotting.md` shows static maps, integration-horizon sliders,
  generated GIFs, and globe plots.

## Basic usage

```julia
using CairoMakie
using Random
using RingGrids
using SpeedyWeatherFTLE

Random.seed!(42)

grid = FullGaussianGrid(8)
u = 25 * rand(grid)
v = 25 * rand(grid)

result = positive_FTLE(
    u,
    v;
    simulation_days = 1,
    dynamics = false,
    rint_hours = 6,
    return_result = true,
    time_indices = :nonzero,
)

fig, ax, sp, cb = surface_plot(result; coastlines = false)
```

This example uses random abstract velocity fields, so the plot disables
coastlines. For geophysical fields, leave coastlines on or style them as visual
context.

Use `negative_FTLE` for backward-time FTLE in frozen prescribed flows. Pass
`time_indices = :last` or `:final` when only the final tracker sample is
needed, or `:nonzero` to skip the initial `0 h` sample where FTLE is undefined.
The high-level prescribed-field API does not yet provide a true unsteady
negative-time FTLE workflow that replays the reversed velocity history of an
evolving SpeedyWeather simulation.

Saved particle files can be post-processed without rerunning the simulation:

```julia
result = get_FTLE(
    u,
    v;
    return_result = true,
    keep_particle_file = true,
    particle_tracker_path = "particle_output",
)

FTLE, time_hours = FTLE_from_particle_file(
    result.particle_file_path,
    result.spectral_grid,
    result.dist_km;
    time_indices = :last,
)
```

## Setting up the project for development
To set up the project for development for the first time, first clone this
repository. Then, from the repository directory, open `julia` and run:

```julia
] activate .
] instantiate
```

`activate .` selects the project in the present directory. `instantiate`
creates or updates the environment needed for local package work; you usually
only need to instantiate once.

## Running the test suite located in `./test`
Run the following code:
```julia
]test
```
