module SpeedyWeatherFTLE

include("visualization/surface_plot.jl")
include("visualization/slider_plot.jl")
include("get_FTLE.jl")

export surface_plot
export slider_plot
export get_FTLE
export Re

end
