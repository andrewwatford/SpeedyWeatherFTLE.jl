# Running FTLE

The main entry point is [`FTLE`](@ref). It accepts zonal and meridional
`RingGrids.Field` values, runs frozen prescribed-flow particle advection
internally, and post-processes the saved deformation into FTLE fields.

```@example simulation
using Logging
using RingGrids
using SpeedyWeatherFTLE

grid = FullGaussianGrid(4)
u = 0 .* rand(grid)
v = 0 .* rand(grid)

ftle, time_hours = with_logger(NullLogger()) do
    FTLE(
        u,
        v;
        simulation_days = 0.25,
        rint_hours = 3,
        particle_advection_every_n_time_steps = 1,
        time_indices = :last,
    )
end

(ftle isa Field, size(ftle), time_hours)
```

`ftle` is a `RingGrids.Field` time series. Use ordinary field indexing to pick
one saved horizon:

```@example simulation
final_field = ftle[:, end]
stretch = stretching_factor(ftle, time_hours)

(final_field isa Field, stretch isa Field)
```

For backward-time FTLE, keep the same function and set `backwards = true`:

```@example simulation
backward_ftle, backward_time_hours = with_logger(NullLogger()) do
    FTLE(
        u,
        v;
        backwards = true,
        simulation_days = 0.25,
        rint_hours = 3,
        particle_advection_every_n_time_steps = 1,
        time_indices = :last,
    )
end

(backward_ftle isa Field, backward_time_hours)
```

`FTLE` uses `radius = SpeedyWeather.DEFAULT_RADIUS` by default. Pass
`radius = ...` in metres to use a different spherical distance scale for both
the internal SpeedyWeather model and the FTLE latitude/longitude conversions.

Dynamic workflows are out of scope for SpeedyWeatherFTLE. The internal
SpeedyWeather simulation used by [`FTLE`](@ref) is always initialized with
`dynamics = false`.
