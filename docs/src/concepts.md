# Concepts and Data Layout

## Finite-Time Lyapunov Exponents

An FTLE measures the largest finite-time stretching rate near each initial
position. SpeedyWeatherFTLE estimates that stretching by releasing four
particles around every grid point: one east, one west, one north, and one south.
Those trajectories define a centred finite-difference approximation to the
displacement gradient. The largest eigenvalue of the resulting right
Cauchy-Green tensor gives the FTLE.

FTLE values returned by this package are in inverse hours because the particle
output times are converted to hours before post-processing.

Use [`stretching_factor`](@ref) to convert FTLE values back to the finite-time
stretching factor `exp(FTLE * T)`. This dimensionless value is often easier to
interpret: `2` means nearby particles separated by a factor of two over the
selected integration time.

## Positive and Negative Direction

Use [`positive_FTLE`](@ref) for forward-time FTLE and [`negative_FTLE`](@ref) for
backward-time FTLE. Both wrappers call [`get_FTLE`](@ref) and fix the
`backwards` keyword for you.

```julia
positive = positive_FTLE(u, v; return_result = true)
negative = negative_FTLE(u, v; return_result = true)
```

When you need to choose the direction programmatically, call [`get_FTLE`](@ref)
directly:

```julia
result = get_FTLE(u, v; backwards = true, return_result = true)
```

## Array Shapes

The main numerical output is an `FTLE_grid_time` matrix with shape
`(grid point, selected time)`. The grid-point order is the order used by
`RingGrids.get_londlatds(grid)`.

For lower-level post-processing, particle longitude and latitude arrays should
have shape `(particle, time)`. There must be four particles per FTLE grid point,
ordered east, west, north, south.

## Selecting Output Times

The `time_indices` keyword controls which particle-tracker time columns are
post-processed. It is available in [`get_FTLE`](@ref),
[`FTLE_from_particle_file`](@ref), and [`FTLE_from_particles`](@ref).

Common selectors are:

- `:nonzero`: skip the initial zero-duration sample where FTLE is undefined.
- `:last` or `:final`: compute only the final tracker sample.
- `:` or `:all`: compute all tracker samples.
- an integer, integer vector/range, or boolean mask: compute explicit samples.

## FTLE Integration Horizons

When you keep multiple tracker output times, the columns of `FTLE_grid_time`
are FTLE estimates for different integration durations from the same particle
release. This is not a conventional time series of instantaneous FTLE fields.
FTLE is a finite-window diagnostic: each column answers what stretching rate is
inferred over an interval of length `T`, for example after integrating for 6 h,
12 h, 18 h, and so on.

The slider plots in these docs therefore sweep integration horizons. They do
not show a single FTLE field evolving in time. A genuine evolving-flow
comparison needs a clear choice of release time or averaging window for each
frame, and is better represented by separate FTLE calculations with those
choices made explicitly.

If you want to compare materially different flow states, run separate FTLE
experiments, for example one initialized from a summer-like jet and one from a
winter-like jet, or one with `dynamics = false` and one with `dynamics = true`.

## Coastlines and Land

FTLE is computed from the supplied velocity field and particle trajectories.
The plotting helpers can draw GeoMakie coastlines, but those coastlines are only
visual overlays. They do not define land, mask output values, or prevent
particles from crossing a shoreline.

If land matters for your experiment, encode it in the model setup, velocity
fields, or post-processing mask. For idealized global examples, coastlines are
best read as orientation marks; for abstract test flows, turn them off with
`coastlines = false`.

## Exact In-Memory Example

For a linear map `A = [2 0; 0 0.5]` over four hours, the largest singular value
is `2`, so the FTLE is `log(2) / 4`.

```@example exact
using SpeedyWeatherFTLE

dist_km = 10
T = 4.0
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
     delta      0.5 * delta
    -delta     -0.5 * delta
]

time_hours = [0.0, T]

ftle, selected_times = FTLE_from_particles(
    plonds_time,
    platds_time,
    time_hours,
    1,
    dist_km;
    time_indices = :nonzero,
)

round(ftle[1, 1], digits = 6), round(log(2) / T, digits = 6)
```

The initial time column is skipped here because FTLE is undefined at `t = 0`.
If you include it, the package writes `NaN` for that sample.
