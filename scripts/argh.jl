using SpeedyWeather

include("lyapunov_functions.jl")

gFTLE, ggrid, time_hours = lyapunov_FTLE(; dynamics=false, backwards=false, use_initial=true, rint_hours=1)

# Example - plot heatmap of FTLE at fixed point in time
FTLE_field = Field(gFTLE[:,end], ggrid) # Select FTLE at end of simulation
#heatmap(FTLE_field)

# Example - plot FTLE as a function of time
# lines(time_hours,gFTLE[2000,:])
