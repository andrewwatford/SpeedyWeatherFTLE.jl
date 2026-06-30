# Running Simulations

The high-level simulation API starts from two velocity fields:

- `u`: zonal velocity on a `RingGrids` grid.
- `v`: meridional velocity on the same grid.

The fields are passed to SpeedyWeather, particles are released around every grid
point, a `ParticleTracker` writes trajectories, and SpeedyWeatherFTLE
post-processes those trajectories into FTLE values.

## Positive-Time FTLE

```julia
using RingGrids
using SpeedyWeatherFTLE

spatial_grid = FullGaussianGrid(20)
u = 100 * rand(spatial_grid)
v = 100 * rand(spatial_grid)

result = positive_FTLE(
    u,
    v;
    simulation_days = 10,
    dynamics = true,
    rint_hours = 3,
    return_result = true,
    time_indices = :nonzero,
)
```

`result.ftle` has one row per grid point and one column per selected output
time. Use [`final_ftle`](@ref) or [`final_ftle_field`](@ref) for the final
selected time.

```julia
ftle_vector = final_ftle(result)
ftle_field_for_plotting = final_ftle_field(result)
```

## Negative-Time FTLE

Negative-time FTLE uses the same velocity fields but runs particle advection
backward in time:

```julia
negative = negative_FTLE(
    u,
    v;
    simulation_days = 10,
    dynamics = true,
    rint_hours = 3,
    return_result = true,
    time_indices = :last,
)
```

Use `time_indices = :last` when you only need the final output. This avoids
post-processing intermediate tracker columns.

## Tuple Return Mode

When `return_result = false`, [`get_FTLE`](@ref) and the direction-specific
wrappers return a tuple:

```julia
FTLE_grid_time, spectral_grid, time_hours = positive_FTLE(
    u,
    v;
    simulation_days = 2,
    rint_hours = 6,
    time_indices = :nonzero,
)
```

This is useful when you want a minimal return value or when existing code
already expects the older tuple interface.

## Direction as a Keyword

Use [`get_FTLE`](@ref) directly when the direction is a runtime choice:

```julia
backward_run = true

result = get_FTLE(
    u,
    v;
    backwards = backward_run,
    return_result = true,
    time_indices = :nonzero,
)
```

The wrappers [`positive_FTLE`](@ref) and [`negative_FTLE`](@ref) intentionally
reject a `backwards` keyword so that the direction cannot be changed by
accident.

## Particle Tracker Controls

The particle file is normally temporary and deleted after FTLE computation. You
can still control the tracker output settings:

```julia
result = positive_FTLE(
    u,
    v;
    return_result = true,
    particle_advection_every_n_time_steps = 6,
    particle_tracker_keepbits = 15,
    particle_tracker_compression_level = 1,
    particle_tracker_shuffle = false,
)
```

Use the [Particle Files](particle_files.md) workflow when you want to keep and
reuse the saved trajectories.
