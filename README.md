# SpeedyWeatherFTLE

[![Build Status](https://github.com/andrewwatford/SpeedyWeatherFTLE.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/andrewwatford/SpeedyWeatherFTLE.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![docs](https://img.shields.io/badge/documentation-latest_release-blue.svg)](https://andrewwatford.github.io/SpeedyWeatherFTLE.jl/)

SpeedyWeatherFTLE computes finite-time Lyapunov exponents (FTLEs) from
SpeedyWeather particle trajectories. It can run a SpeedyWeather particle
tracking simulation from prescribed velocity fields, compute positive- or
negative-time FTLE, reuse saved `ParticleTracker` NetCDF files, and convert
FTLE arrays to RingGrids fields for plotting.

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
julia --project=docs docs/make.jl
```

## Examples

See `examples/speedyweather_ftle_snapshots.ipynb` for a local notebook that runs
and caches a dynamic SpeedyWeather flow simulation, plots the
initial and evolved velocity components, and then compares initial and evolved
FTLE fields for a user-selected integration horizon.

## Basic usage

```julia
using SpeedyWeatherFTLE, RingGrids

grid = FullGaussianGrid(20)
u = 100 * rand(grid)
v = 100 * rand(grid)

result = positive_FTLE(
    u,
    v;
    dynamics = true,
    return_result = true,
    time_indices = :nonzero,
)

field = final_ftle_field(result)
fig, ax, sp, cb = surface_plot(final_ftle(result), result.spectral_grid)
```

Use `negative_FTLE` for backward-time FTLE. Pass `time_indices = :last` or
`:final` when only the final tracker sample is needed, or `:nonzero` to skip
the initial `0 h` sample where FTLE is undefined.

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
To set up the project for development for the first time, first clone this repository. Then, from the repository directory, open `julia` and run:
```julia
]instantiate
]activate .
```
`instantiate` creates a `Manifest.toml` file, which creates the environment that is needed for all future use of the package - you usually only need to instantiate once. `activate .` activates the environment in the present directory, that is, the current project.

## Running the test suite located in `./test`
Run the following code:
```julia
]test
```
