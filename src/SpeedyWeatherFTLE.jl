module SpeedyWeatherFTLE

include("mock/greet.jl")
include("visualization/surface_plot.jl")
include("lyapunov_FTLE.jl")

export greet
export surface_plot
export lyapunov_FTLE
export Re

end
