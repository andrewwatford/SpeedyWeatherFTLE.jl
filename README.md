# SpeedyWeatherFTLE

SpeedyWeatherFTLE computes finite-time Lyapunov exponents (FTLEs) from
SpeedyWeather particle trajectories. The package is intentionally small:

- build the four-particle east/west/north/south release stencil;
- compute FTLE from in-memory particle positions;
- compute FTLE from `ParticleTracker` NetCDF files;
- optionally run a frozen prescribed-velocity SpeedyWeather experiment;
- convert FTLE arrays to `RingGrids.Field` for your own plotting code.

There is no bundled plotting layer.

## Installation

Requires Julia 1.12 or newer.

SpeedyWeatherFTLE currently needs particle-tracking APIs from the
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

## Particle Trajectories

`FTLE_from_particles` is the smallest workflow. `plonds_time` and `platds_time`
must be shaped `(particle, time)`, with four particles per grid point in east,
west, north, south order.

```julia
using SpeedyWeatherFTLE

dist_km = 10.0
delta = rad2deg(dist_km * 1000 / SpeedyWeatherFTLE.Re)

plonds_time = [
     delta      2 * delta
    -delta     -2 * delta
     0.0        0.0
     0.0        0.0
]

platds_time = [
     0.0        0.0
     0.0        0.0
     delta      delta
    -delta     -delta
]

ftle, time_hours = FTLE_from_particles(
    plonds_time,
    platds_time,
    [1.0, 2.0],
    1,
    dist_km;
    time_indices = :last,
)
```

Use `initial_FTLE_particle_positions(grid, dist_km)` to create the matching
release stencil for a `RingGrids` grid or `SpeedyWeather.SpectralGrid`.

## Particle Files

Saved SpeedyWeather particle files can be post-processed without rerunning the
simulation:

```julia
ftle, time_hours = FTLE_from_particle_file(
    "particles.nc",
    spectral_grid,
    10.0;
    time_indices = :nonzero,
)
```

The file must contain `time`, `lon`, and `lat` variables with longitude and
latitude dimensions `(particle, time)`.

## Frozen Fields

For a prescribed frozen velocity field, call `positive_FTLE` or
`negative_FTLE`:

```julia
using RingGrids
using SpeedyWeatherFTLE

grid = FullGaussianGrid(8)
u = rand(grid)
v = rand(grid)

result = positive_FTLE(
    u,
    v;
    simulation_days = 1,
    rint_hours = 6,
    return_result = true,
    time_indices = :nonzero,
)

final = final_ftle(result)
field = final_ftle_field(result)
stretch = stretching_factor(result)
```

`get_FTLE(u, v; backwards = true)` is the lower-level form. The wrapper supports
frozen prescribed fields only; pass `dynamics = false`.

## Public API

- `initial_FTLE_particle_positions`, `initial_FTLE_particle_positions!`
- `FTLE_from_particles`, `FTLE_from_particles!`
- `FTLE_from_particle_file`, `FTLE_from_particle_file!`
- `get_FTLE`, `positive_FTLE`, `negative_FTLE`
- `FTLEResult`, `final_ftle`, `final_ftle_field`, `ftle_field`
- `stretching_factor`, `stretching_factor!`

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
