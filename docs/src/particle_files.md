# Particle Files

SpeedyWeatherFTLE computes FTLE from the same NetCDF files written by
SpeedyWeather's `ParticleTracker`. Saving those files is useful when particle
tracking is expensive and you want to try different post-processing choices
without rerunning the simulation.

## Keep the File from a Simulation

Set `keep_particle_file = true` and `return_result = true` to store the file
path in the returned [`FTLEResult`](@ref):

```julia
using RingGrids
using SpeedyWeatherFTLE

spatial_grid = FullGaussianGrid(20)
u = 100 * rand(spatial_grid)
v = 100 * rand(spatial_grid)

result = positive_FTLE(
    u,
    v;
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
displacement-gradient work array once, then call
[`FTLE_from_particle_file!`](@ref):

```julia
npoints = result.spectral_grid.npoints
FTLE_buffer = Matrix{Float64}(undef, npoints, length(result.time_hours))
B = Array{Float64}(undef, 2, 2, npoints)

FTLE_from_particle_file!(
    FTLE_buffer,
    B,
    result.particle_file_path,
    result.spectral_grid,
    result.dist_km;
    time_indices = eachindex(result.time_hours),
)
```

When `time_indices` selects fewer columns, size `FTLE_buffer` with one column per
selected time.

## In-Memory Trajectories

If you already have particle trajectories in memory, call
[`FTLE_from_particles`](@ref) directly. Longitude and latitude arrays must be
shaped as `(particle, time)`, with four particles per FTLE grid point.

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
