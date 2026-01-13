# Plotting the Finite-Time Lyapunov Exponent (FTLE) for climatological initial conditions

In this example we plot the forwards-time finite-time Lyapunov exponents (FTLEs) using the default climatological initial conditions in SpeedyWeather. This is the default behaviour of our `lyapunov_FTLE` function.

```@example
using SpeedyWeatherFTLE
pFTLE, grid, time_hours = lyapunov_FTLE()
pFTLE_final = Field(pFTLE[:, end], grid)
fig, ax, sp, cb = surface_plot(
                pFTLE_final;
                title = "pFTLE field at $(time_hours[end]) hours",
                label = "pFTLE [1/s]",
                )
fig
```
