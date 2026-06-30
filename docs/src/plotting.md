# Plotting

The plotting helpers convert FTLE vectors or time-series matrices into
`RingGrids.Field` objects and then interpolate them onto a regular longitude
and latitude grid for GeoMakie.

Use CairoMakie in scripts and documentation builds, or GLMakie locally when you
want interactive windows.

## Plot One Output Time

[`surface_plot`](@ref) accepts an [`FTLEResult`](@ref), an FTLE vector plus a
grid, an FTLE matrix plus a grid, or a `RingGrids.Field`.

```@example plotting
using CairoMakie
using RingGrids
using SpeedyWeatherFTLE

spatial_grid = FullGaussianGrid(8)
londs, latds = RingGrids.get_londlatds(spatial_grid)

synthetic_ftle = @. 0.015 + 0.005 * sind(latds)^2

fig, ax, sp, cb = surface_plot(
    synthetic_ftle,
    spatial_grid;
    title = "Synthetic FTLE",
    label = "FTLE [1/h]",
)

fig
```

If you have an [`FTLEResult`](@ref), this is enough:

```julia
fig, ax, sp, cb = surface_plot(result; label = "FTLE [1/h]")
```

## Plot a Time Series

[`slider_plot`](@ref) accepts the result object directly or the tuple-style
output from [`get_FTLE`](@ref).

```@example plotting
time_hours = [6.0, 12.0, 18.0]
FTLE_grid_time = hcat(synthetic_ftle, 1.5 .* synthetic_ftle, 2 .* synthetic_ftle)

fig, ax, sp, cb = slider_plot(
    time_hours,
    FTLE_grid_time,
    spatial_grid;
    title = "Synthetic FTLE time series",
    colorbar_label = "FTLE [1/h]",
)

fig
```

For result objects:

```julia
fig, ax, sp, cb = slider_plot(
    result;
    title = "FTLE time series",
    colorbar_label = "FTLE [1/h]",
)
```

The zero-duration tracker sample is skipped by default because FTLE is
undefined at `t = 0`. Pass `start_index = 1` if you explicitly want to include
that column.

## Convert Without Plotting

Use [`ftle_field`](@ref) when you want a `RingGrids.Field` for your own plotting
or interpolation code:

```julia
field_ts = ftle_field(FTLE_grid_time, spatial_grid)
field_final = ftle_field(FTLE_grid_time[:, end], spatial_grid)
field_from_result = ftle_field(result; time_indices = :last)
```
