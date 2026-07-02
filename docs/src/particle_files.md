# Particle Files

SpeedyWeatherFTLE computes FTLE from FTLE-compatible NetCDF files written by
SpeedyWeather's `ParticleTracker`. Saving those files is useful when particle
tracking is expensive and you want to try different post-processing choices
without rerunning the simulation.

Not every `ParticleTracker` file is FTLE-compatible. Before post-processing a
file, make sure it satisfies the same release layout used by
[`initial_FTLE_particle_positions`](@ref):

- four particles per FTLE grid point;
- particle order east, west, north, south for every grid point;
- initial positions generated on the same grid or `SpectralGrid` passed to
  [`FTLE_from_particle_file`](@ref);
- the same `dist_km` used for release and post-processing;
- elapsed output times that represent integration duration from the release.

Files created by [`positive_FTLE`](@ref), [`negative_FTLE`](@ref), or
[`get_FTLE`](@ref) with `keep_particle_file = true` use this stencil and store
SpeedyWeatherFTLE metadata for the requested `dist_km`, particle order, and
particle-tracker keepbits. For files from custom SpeedyWeather simulations,
prefer generating the initial particle positions with
[`initial_FTLE_particle_positions`](@ref) so the file layout and the
post-processing assumptions stay matched.

For an existing SpeedyWeather `Simulation`, use [`attach_FTLE_tracker!`](@ref)
instead of wiring the particle positions and tracker by hand. It prepares the
FTLE particle stencil, attaches a `ParticleTracker`, and returns an
[`FTLEParticleSetup`](@ref) containing the `spectral_grid`, `dist_km`, particle
order, output path, and elapsed-time semantics used by
[`FTLE_from_particle_file`](@ref).

By default, [`FTLE_from_particle_file`](@ref) and
[`FTLE_from_particle_file!`](@ref) validate that the file has the required
`particle` and `time` dimensions, `time`, `lon`, and `lat` variables, position
variables shaped as `(particle, time)`, the expected particle count, matching
SpeedyWeatherFTLE metadata when present, and an initial longitude/latitude
column that matches the canonical east, west, north, south stencil for the
supplied grid or `SpectralGrid` and `dist_km`. This catches files produced with
a different grid, stencil order, or perturbation distance before FTLE is
computed.

The initial-position tolerance is intentionally small enough that a default
`dist_km = 10` file is rejected if post-processed as `dist_km = 20`. Very
coarsely quantized legacy files may fail strict validation because the saved
coordinates no longer distinguish nearby stencils. Regenerate those files with
SpeedyWeatherFTLE metadata when possible; otherwise disable initial-position
validation only after checking compatibility externally.

Only opt out when you have already validated compatibility or intentionally
need structural post-processing without the initial-position check:

```julia
FTLE_grid_time, time_hours = FTLE_from_particle_file(
    path,
    spectral_grid,
    dist_km;
    time_indices = :nonzero,
    validate_initial_positions = false,
)
```

With validation disabled, the file must still have compatible `lon`, `lat`, and
`time` variables and four particles per FTLE grid point. Non-finite particle
output times or selected particle coordinates throw an `ArgumentError`;
zero-duration samples are returned as `NaN` if selected.

## Keep the File from a Simulation

Set `keep_particle_file = true` and `return_result = true` to store the file
path in the returned [`FTLEResult`](@ref):

```julia
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
    keep_particle_file = true,
    particle_tracker_path = "particle_output",
    particle_tracker_filename = "particles.nc",
    time_indices = :last,
)

result.particle_file_path
```

The file is deleted automatically unless `keep_particle_file = true` or
`return_particle_file_path = true`.

## Reprocess a Saved File

Use [`FTLE_from_particle_file`](@ref) to compute FTLE from an existing
`ParticleTracker` file:

```julia
FTLE_grid_time, time_hours = FTLE_from_particle_file(
    result.particle_file_path,
    result.spectral_grid,
    result.dist_km;
    time_indices = :nonzero,
)
```

When the file was created with [`attach_FTLE_tracker!`](@ref), pass the setup
directly:

```julia
FTLE_grid_time, time_hours = FTLE_from_particle_file(setup)
```

You can also post-process only the final saved time:

```julia
FTLE_final, final_time_hours = FTLE_from_particle_file(
    result.particle_file_path,
    result.spectral_grid,
    result.dist_km;
    time_indices = :last,
)
```

## Reuse Allocations

For repeated post-processing, allocate the output matrix and the
deformation-gradient work array once, then call
[`FTLE_from_particle_file!`](@ref):

```julia
npoints = result.spectral_grid.npoints
FTLE_buffer = Matrix{Float64}(undef, npoints, 1)
B = Array{Float64}(undef, 2, 2, npoints)

FTLE_from_particle_file!(
    FTLE_buffer,
    B,
    result.particle_file_path,
    result.spectral_grid,
    result.dist_km;
    time_indices = :last,
)
```

When `time_indices` selects fewer columns, size `FTLE_buffer` with one column per
selected time. `validate_initial_positions = true` is also the default for the
in-place form; pass `validate_initial_positions = false` only for a
pre-validated compatible file.

## In-Memory Trajectories

If you already have particle trajectories in memory, call
[`FTLE_from_particles`](@ref) directly. Longitude and latitude arrays must be
shaped as `(particle, time)`, with four particles per FTLE grid point.
Sampling times are elapsed durations in hours; signed values are accepted and
the FTLE rate uses their magnitude, while non-finite times throw.

Use [`initial_FTLE_particle_positions`](@ref) to create the initial east, west,
north, south release stencil for your own trajectory integrator:

```julia
initial_lon, initial_lat = initial_FTLE_particle_positions(spectral_grid, dist_km)
```

```julia
FTLE_grid_time, selected_time_hours = FTLE_from_particles(
    particle_lon,
    particle_lat,
    time_hours,
    spectral_grid,
    dist_km;
    time_indices = :nonzero,
)
```

Use [`FTLE_from_particles!`](@ref) when you want to provide the output and work
arrays yourself.
