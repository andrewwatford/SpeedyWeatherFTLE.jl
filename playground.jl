using SpeedyWeather, GeoMakie, CairoMakie

# Set up a random u and v field on a FullGaussianGrid
u_data = 1000 * randn(Float64, 8, 4)
u_field = FullGaussianGrid(u_data, input_as=Matrix)
v_data = 1000 * randn(Float64, 8, 4)
v_field = FullGaussianGrid(v_data, input_as=Matrix)

# Run FTLE calculation
pFTLE, grid, time_hours = lyapunov_FTLE(
    trunc=35,
    use_climatological=false,
    use_random=false,
    zonal_velocity_field=u_field,
    meridional_velocity_field=v_field,
    rint_hours=12,
)
pFTLE_final = Field(pFTLE[:, end], grid)
fig, ax, sp, cb = surface_plot(
                pFTLE_final;
                title = "pFTLE field at $(time_hours[end]) hours",
                label = "pFTLE [1/h]",
                )
fig