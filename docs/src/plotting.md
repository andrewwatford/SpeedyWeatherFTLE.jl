# Plotting

The plotting helpers are optional. They load through a package extension after
GeoMakie is available, so compute-only users can load SpeedyWeatherFTLE without
loading Makie.

```@example plotting
using CairoMakie
using GeoMakie
using RingGrids
using SpeedyWeather
using SpeedyWeatherFTLE

CairoMakie.activate!()

spectral_grid = SpectralGrid(nlayers = 1, trunc = 8, Grid = FullGaussianGrid)
values = zeros(spectral_grid.npoints, 3)
values[:, 1] .= NaN
values[:, 2] .= range(0.0, 0.20; length = spectral_grid.npoints)
values[:, 3] .= range(0.05, 0.30; length = spectral_grid.npoints)

ftle = Field(values, spectral_grid.grid)
time_hours = [0.0, 6.0, 12.0]

nothing # hide
```

In normal use, `ftle` and `time_hours` come from [`FTLE`](@ref). The synthetic
field above keeps the plotting examples small and reproducible.

## Surface Plot

[`surface_plot`](@ref) accepts a `RingGrids.Field` directly. For a time series
field, it uses the final selected column by default and labels the colorbar as
`FTLE [1/h]`.

```@example plotting
fig, ax, sp, cb = surface_plot(
    ftle;
    time_hours,
    coastlines = false,
    title = "Forward FTLE",
)

nothing # hide
```

Select a saved integration horizon by index or by nearest hour:

```@example plotting
surface_plot(ftle; time_hours, time_index = 3, coastlines = false);
surface_plot(ftle; time_hours, time_hour = 12, coastlines = false);

nothing # hide
```

Use [`shared_colorrange`](@ref) for comparable plots:

```@example plotting
shared = shared_colorrange(ftle; pad = 0.05)
surface_plot(ftle; time_hours, colorrange = shared, coastlines = false);

nothing # hide
```

## Slider Plot

[`slider_plot`](@ref) shows each saved nonzero integration horizon. It skips
zero-duration FTLE columns by default because FTLE is undefined at release.

```@example plotting
handle = slider_plot(
    time_hours,
    ftle;
    return_handle = true,
    coastlines = false,
)

set_slider_time!(handle, 12)

nothing # hide
```

The returned [`SliderPlotHandle`](@ref) exposes the figure, axis, surface,
colorbar, slider, live label, and plotted times.

Record an animation with the same controls:

```@example plotting
fake_record(callback, fig, path, frames; framerate, kwargs...) = (foreach(callback, frames); path) # hide

animate_slider_plot(
    "ftle.gif",
    time_hours,
    ftle;
    framerate = 6,
    coastlines = false,
    record_function = fake_record, # hide
)
```

## Globe Plot

[`globe_plot`](@ref) uses GeoMakie's `GlobeAxis`. With `GLMakie` active, the
globe can be rotated interactively. With `CairoMakie`, it renders as a static
figure.

```@example plotting
fig, ax, sp, cb = globe_plot(
    ftle;
    time_hours,
    time_hour = 12,
    coastlines = false,
    colorbar = true,
)

nothing # hide
```
