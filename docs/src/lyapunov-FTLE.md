# Plotting the Finite-Time Lyapunov Exponent (FTLE)

In this example we plot the forwards-time finite-time Lyapunov exponents (FTLEs) using the default climatological initial conditions in SpeedyWeather. This is the default behaviour of our `lyapunov_FTLE` function.

```@example
using SpeedyWeatherFTLE
using RingGrids: Field
pFTLE, grid, time_hours = lyapunov_FTLE()
pFTLE_final = Field(pFTLE[:, end], grid)
fig, ax, sp, cb = surface_plot(
                pFTLE_final;
                title = "pFTLE field at $(time_hours[end]) hours",
                label = "pFTLE [1/h]",
                )
fig
```

We can also use random initial conditions:

```@example
using SpeedyWeatherFTLE
using RingGrids: Field
pFTLE, grid, time_hours = lyapunov_FTLE(use_random=true, use_climatological=false)
pFTLE_final = Field(pFTLE[:, end], grid)
fig, ax, sp, cb = surface_plot(
                pFTLE_final;
                title = "pFTLE field at $(time_hours[end]) hours",
                label = "pFTLE [1/h]",
                )
fig
```

And, finally, we can prescribe our own:
```@example
using SpeedyWeatherFTLE
using RingGrids: Field, FullGaussianGrid

# Set up a random u and v field on a FullGaussianGrid
u_data = 1000 * randn(Float64, 8, 4)
u_field = FullGaussianGrid(u_data, input_as=Matrix)
v_data = 1000 * randn(Float64, 8, 4)
v_field = FullGaussianGrid(v_data, input_as=Matrix)

pFTLE, grid, time_hours = lyapunov_FTLE(use_random=false, use_climatological=false, zonal_velocity_field=u_field, meridional_velocity_field=v_field)
pFTLE_final = Field(pFTLE[:, end], grid)
fig, ax, sp, cb = surface_plot(
                pFTLE_final;
                title = "pFTLE field at $(time_hours[end]) hours",
                label = "pFTLE [1/h]",
                )
fig
```
