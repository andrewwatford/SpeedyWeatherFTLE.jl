# Plotting

The plotting helpers are optional. They load through a package extension after
GeoMakie is available, so compute-only users can load SpeedyWeatherFTLE without
loading Makie.

```julia
using CairoMakie
using GeoMakie
using RingGrids
using SpeedyWeather
using SpeedyWeatherFTLE

CairoMakie.activate!()

spectral_grid = SpectralGrid(nlayers = 1, trunc = 8, Grid = FullGaussianGrid)
ftle = zeros(spectral_grid.npoints, 3)
ftle[:, 1] .= NaN
ftle[:, 2] .= range(0.0, 0.20; length = spectral_grid.npoints)
ftle[:, 3] .= range(0.05, 0.30; length = spectral_grid.npoints)

result = FTLEResult(
    ftle,
    spectral_grid,
    [0.0, 6.0, 12.0];
    dist_km = 10,
    backwards = false,
    rint_hours = 6,
)
```

In normal use, `result` comes from [`FTLE`](@ref) with `return_result = true`.
The synthetic result above keeps the plotting examples small and reproducible.

## Surface Plot

[`surface_plot`](@ref) accepts an [`FTLEResult`](@ref) directly. By default it
uses the final selected FTLE column, labels the colorbar as `FTLE [1/h]`, and
uses finite FTLE values for the color range.

```julia
fig, ax, sp, cb = surface_plot(
    result;
    coastlines = false,
    title = "Forward FTLE",
)
```

Select a saved integration horizon by index or by nearest hour:

```julia
surface_plot(result; time_index = 3, coastlines = false)
surface_plot(result; time_hour = 12, coastlines = false)
```

The same helper accepts raw FTLE arrays when you pass the grid:

```julia
surface_plot(result.ftle, result.spectral_grid; time_hours = result.time_hours, time_hour = 12)
surface_plot(final_ftle(result), result.spectral_grid)
```

Use [`ftle_colorrange`](@ref) for comparable plots:

```julia
shared = ftle_colorrange(result; pad = 0.05)
backward = FTLEResult(
    reverse(result.ftle; dims = 1),
    result.spectral_grid,
    result.time_hours;
    dist_km = result.dist_km,
    backwards = true,
    rint_hours = result.rint_hours,
)

surface_plot(result; colorrange = shared, coastlines = false)
surface_plot(backward; colorrange = shared, coastlines = false)
```

## Slider Plot

[`slider_plot`](@ref) shows each saved nonzero integration horizon. It skips
zero-duration FTLE columns by default because FTLE is undefined at release.

```julia
handle = slider_plot(
    result;
    return_handle = true,
    coastlines = false,
)

set_slider_time!(handle, 12)
```

The returned [`SliderPlotHandle`](@ref) exposes the figure, axis, surface,
colorbar, slider, live label, and plotted times.

Record an animation with the same controls:

```julia
animate_slider_plot(
    "ftle.gif",
    result;
    framerate = 6,
    coastlines = false,
)
```

## Globe Plot

[`globe_plot`](@ref) uses GeoMakie's `GlobeAxis`. With `GLMakie` active, the
globe can be rotated interactively. With `CairoMakie`, it renders as a static
figure.

```julia
fig, ax, sp, cb = globe_plot(
    result;
    time_hour = 12,
    coastlines = false,
    colorbar = true,
)
```

`surface_plot` and `globe_plot` accept `RingGrids.Field`, FTLE vectors, FTLE
matrices with a grid, and [`FTLEResult`](@ref). `slider_plot` and
`animate_slider_plot` accept time vectors with time-dependent fields, FTLE
matrices with a grid, and [`FTLEResult`](@ref).
