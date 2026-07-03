# Running FTLE

The main entry point is [`FTLE`](@ref). It accepts zonal and meridional
`RingGrids.Field` values, runs frozen prescribed-flow particle advection
internally, and post-processes the saved deformation into FTLE values.

```julia
using Random
using RingGrids
using SpeedyWeatherFTLE

Random.seed!(42)

grid = FullGaussianGrid(8)
u = 25 .* rand(grid)
v = 25 .* rand(grid)

result = FTLE(
    u,
    v;
    simulation_days = 1,
    rint_hours = 6,
    particle_advection_every_n_time_steps = 1,
    return_result = true,
    time_indices = :nonzero,
)
```

`result` is an [`FTLEResult`](@ref):

```julia
result.ftle          # matrix with dimensions (grid point, selected time)
result.spectral_grid # SpeedyWeather SpectralGrid used for the run
result.time_hours    # selected integration horizons
result.direction     # :forward or :backward
```

Common follow-up operations:

```julia
final_values = final_ftle(result)
final_field = final_ftle_field(result)
field_at_12h = ftle_field(result; time_hour = 12)
stretch = stretching_factor(result)
```

For backward-time FTLE, keep the same function and set `backwards = true`:

```julia
backward = FTLE(
    u,
    v;
    backwards = true,
    simulation_days = 1,
    rint_hours = 6,
    particle_advection_every_n_time_steps = 1,
    return_result = true,
    time_indices = :last,
)
```

By default, [`FTLE`](@ref) returns `(ftle, spectral_grid, time_hours)` instead
of an [`FTLEResult`](@ref):

```julia
ftle, spectral_grid, time_hours = FTLE(
    u,
    v;
    simulation_days = 1,
    rint_hours = 6,
    time_indices = :last,
)
```

This wrapper is intentionally for frozen supplied velocity fields. Pass
`dynamics = false` or leave it at the default. Dynamic model initialization
belongs in a dedicated SpeedyWeather simulation workflow.
