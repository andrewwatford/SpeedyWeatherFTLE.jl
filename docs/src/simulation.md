# Running Simulations

The high-level simulation API starts from two velocity fields:

- `u`: zonal velocity on a `RingGrids` grid.
- `v`: meridional velocity on the same grid.

The fields are passed to SpeedyWeather, particles are released around every grid
point, a `ParticleTracker` writes trajectories, and SpeedyWeatherFTLE
post-processes those trajectories into FTLE values.

Simulation progress and warnings use SpeedyWeather's normal output behavior, so
longer particle-tracking runs remain transparent in scripts and notebooks.

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

## Comparing Flow States

For a more geophysical workflow, run separate FTLE experiments for separate
flow states. The example below builds two idealized meandering jets, one
weaker and lower-latitude "summer-like" jet and one stronger, sharper
"winter-like" jet. It then compares their final positive-time FTLE after the
same integration duration.

This is a recipe rather than a doctest because it launches full
SpeedyWeather particle-tracking runs. Start with a small grid while tuning
parameters, then increase the resolution and integration time once the workflow
looks right.

First define and inspect the two prescribed flow states. These lightweight
plots are rendered in the docs so you can check that the summer-like and
winter-like jets are meaningfully different before launching particle tracking.

```@example meandering_jets
using CairoMakie
using RingGrids
using SpeedyWeatherFTLE

function meandering_jet(grid; jet_lat, jet_speed, jet_width = 12, wave_number = 5)
    londs, latds = RingGrids.get_londlatds(grid)
    envelope = @. exp(-((latds - jet_lat) / jet_width)^2)
    phase = @. wave_number * londs

    u = rand(grid)
    v = rand(grid)
    u .= @. jet_speed * envelope * (1 + 0.25 * cosd(phase))
    v .= @. 0.20 * jet_speed * envelope * sind(phase)

    return u, v
end

grid = FullGaussianGrid(32)

summer_u, summer_v = meandering_jet(
    grid;
    jet_lat = 35,
    jet_speed = 35,
    jet_width = 14,
)

winter_u, winter_v = meandering_jet(
    grid;
    jet_lat = 45,
    jet_speed = 55,
    jet_width = 10,
)

summer_speed = sqrt.(summer_u.^2 .+ summer_v.^2)
winter_speed = sqrt.(winter_u.^2 .+ winter_v.^2)
speed_colorrange = (0, 70)

fig_summer_speed, ax_summer_speed, sp_summer_speed, cb_summer_speed = surface_plot(
    summer_speed;
    title = "Summer-like meandering jet speed",
    label = "speed [m/s]",
    colormap = :batlow,
    colorrange = speed_colorrange,
    coastline_linewidth = 0.7,
)

fig_summer_speed
```

```@example meandering_jets
fig_winter_speed, ax_winter_speed, sp_winter_speed, cb_winter_speed = surface_plot(
    winter_speed;
    title = "Winter-like meandering jet speed",
    label = "speed [m/s]",
    colormap = :batlow,
    colorrange = speed_colorrange,
    coastline_linewidth = 0.7,
)

fig_winter_speed
```

Then run one FTLE experiment per flow state and compare the final integration
horizon with a shared color scale:

```julia
common_kwargs = (
    simulation_days = 7,
    dist_km = 25,
    rint_hours = 12,
    particle_advection_every_n_time_steps = 1,
    particle_tracker_keepbits = 20,
    return_result = true,
    time_indices = :last,
)

summer = positive_FTLE(summer_u, summer_v; common_kwargs..., dynamics = true)
winter = positive_FTLE(winter_u, winter_v; common_kwargs..., dynamics = true)

shared_colorrange = ftle_colorrange(final_ftle(summer), final_ftle(winter))

fig_summer, ax_summer, sp_summer, cb_summer = surface_plot(
    summer;
    title = "Summer-like jet FTLE after $(only(summer.time_hours)) h",
    colorrange = shared_colorrange,
)

fig_winter, ax_winter, sp_winter, cb_winter = surface_plot(
    winter;
    title = "Winter-like jet FTLE after $(only(winter.time_hours)) h",
    colorrange = shared_colorrange,
)

display(fig_summer)
display(fig_winter)
```

To isolate the effect of SpeedyWeather's evolving dynamics, keep the same
initial velocity field and compare a frozen-flow run against an evolving-flow
run:

```julia
frozen = positive_FTLE(
    winter_u,
    winter_v;
    common_kwargs...,
    dynamics = false,
)

evolving = positive_FTLE(
    winter_u,
    winter_v;
    common_kwargs...,
    dynamics = true,
)

evolution_colorrange = ftle_colorrange(frozen, evolving)

fig, ax, sp, cb = surface_plot(
    evolving;
    title = "Evolving-flow FTLE after $(only(evolving.time_hours)) h",
    colorrange = evolution_colorrange,
)
```

For an integration-duration sweep from a single release, set
`time_indices = :nonzero` and pass the returned result to [`slider_plot`](@ref).
Read that slider as different FTLE integration horizons, not as an independent
time series of flow states.

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
