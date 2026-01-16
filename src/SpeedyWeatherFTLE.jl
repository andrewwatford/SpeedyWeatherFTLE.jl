module SpeedyWeatherFTLE

using SpeedyWeather, RingGrids
using GeoMakie, Makie
using LinearAlgebra
using NCDatasets

const Re = 6.371e6 # Average Earth radius in meters

include("./FTLE_computations.jl")
include("visualization/surface_plot.jl")
include("visualization/slider_plot.jl")
include("get_FTLE.jl")

export surface_plot
export slider_plot
export get_FTLE
export Re

end
