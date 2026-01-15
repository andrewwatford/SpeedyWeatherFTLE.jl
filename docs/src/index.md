# SpeedyWeatherFTLE

Documentation for `SpeedyWeatherFTLE`.

## Plotting the Finite-Time Lyapunov Exponent (FTLE)

Here we show how to obtain the positive- and negative-time Lyapunov exponents (FTLEs) using SpeedyWeatherFTLE's get_FTLE function, and to plot them using the surface_plot function.

```@example ftle
using SpeedyWeatherFTLE, RingGrids

nlat_half = 20
spatial_grid = FullGaussianGrid(nlat_half)
u = 100 * rand(spatial_grid)
v = 100 * rand(spatial_grid)

pFTLE, p_spectral_grid, p_time_hours = get_FTLE(u, v; dynamics=true, backwards=false)
pFTLE_final = Field(pFTLE[:, end], spatial_grid)

nFTLE, n_spectral_grid, n_time_hours = get_FTLE(u, v; dynamics=true, backwards=true)
nFTLE_final = Field(nFTLE[:, end], spatial_grid)

p_fig, p_ax, p_sp, p_cb = surface_plot(
    pFTLE_final;
    title = "pFTLE field at $(p_time_hours[end]) hours",
    label = "pFTLE [1/h]",
)
p_fig
```

```@example ftle
n_fig, n_ax, n_sp, n_cb = surface_plot(
    nFTLE_final;
    title = "nFTLE field at -$(n_time_hours[end]) hours",
    label = "nFTLE [1/h]",
)
n_fig
```