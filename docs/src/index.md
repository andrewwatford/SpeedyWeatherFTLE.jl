# SpeedyWeatherFTLE

Documentation for `SpeedyWeatherFTLE`.

## Plotting the Finite-Time Lyapunov Exponent (FTLE)

Here we show how to obtain the positive- and negative-time Lyapunov exponents (FTLEs) using SpeedyWeatherFTLE's get_FTLE function and to plot them using the surface_plot and slider_plot functions.

```@example ftle
using CairoMakie
using SpeedyWeatherFTLE, RingGrids

nlat_half = 20
spatial_grid = FullGaussianGrid(nlat_half)
u = 100 * rand(spatial_grid)
v = 100 * rand(spatial_grid)

# using dynamics = true evolves the velocity fields, backwards = false for positive-time FTLE
pFTLE, p_spectral_grid, p_time_hours = get_FTLE(u, v; dynamics=true, backwards=false) 
pFTLE_field = Field(pFTLE, spatial_grid)
pFTLE_final = pFTLE_field[:,end]

p_fig, p_ax, p_sp, p_cb = surface_plot(
    pFTLE_final;
    title = "pFTLE field at $(p_time_hours[end]) hours",
    label = "pFTLE [1/h]",
)

p_fig
```

We can also create a slider plot of the whole time-series. (The slider is not interactive in the docs, but using GLMakie locally it is interactive.)

```@example ftle
fig, ax, sp, cb = slider_plot(
                p_time_hours[10:end],
                pFTLE_field[:,10:end]; # remove noisy initial data
                title = "pFTLE field",
                colorbar_label = "pFTLE [1/h]",
            )

fig
```