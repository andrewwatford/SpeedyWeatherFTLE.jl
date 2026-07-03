# SpeedyWeatherFTLE

SpeedyWeatherFTLE computes finite-time Lyapunov exponents (FTLEs) from
SpeedyWeather particle trajectories. It deliberately stays narrow: create the
FTLE particle stencil, post-process tracked particles, and return arrays or
`RingGrids.Field` values that downstream analysis and plotting code can use.

## Install

SpeedyWeatherFTLE requires Julia 1.12 and the `mk/lyapunov2` SpeedyWeather
source branch because that branch provides the particle advection API used by
the package.

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

## Particle Layout

For each FTLE grid point, release four particles in east, west, north, south
order. [`initial_FTLE_particle_positions`](@ref) creates this layout from
longitude/latitude vectors, a `RingGrids` grid, or a `SpeedyWeather.SpectralGrid`.
The longitude and latitude trajectory arrays consumed by
[`FTLE_from_particles`](@ref) must be shaped `(particle, time)`.

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

## Particle Files

Use [`FTLE_from_particle_file`](@ref) for saved SpeedyWeather `ParticleTracker`
NetCDF files. Files must contain `time`, `lon`, and `lat` variables, with
longitude and latitude dimensions `(particle, time)`.

```julia
ftle, time_hours = FTLE_from_particle_file(
    "particles.nc",
    spectral_grid,
    10.0;
    time_indices = :nonzero,
)
```

By default, particle files are checked against the expected release stencil.
Use `validate_initial_positions = false` only for externally validated legacy
files or synthetic tests.

## Frozen Velocity Fields

For a frozen prescribed velocity field, use [`positive_FTLE`](@ref),
[`negative_FTLE`](@ref), or [`get_FTLE`](@ref).

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

final_ftle(result)
final_ftle_field(result)
stretching_factor(result)
```

`get_FTLE` intentionally supports frozen supplied velocity fields only. Dynamic
model initialization belongs in a dedicated SpeedyWeather workflow.

## Development

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
Pkg.test()
```

Build these docs with:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```
