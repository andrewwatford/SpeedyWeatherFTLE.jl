module SpeedyWeatherFTLE

using SpeedyWeather, RingGrids
using GeoMakie, Makie
using LinearAlgebra
using NCDatasets
using Logging

"""
    Re

Average Earth radius used by SpeedyWeatherFTLE, in metres.
"""
const Re = 6.371e6 # Average Earth radius in meters

include("grid_helpers.jl")
include("./FTLE_computations.jl")
include("visualization/ftle_field.jl")
include("FTLE_result.jl")
include("FTLE_diagnostics.jl")
include("visualization/surface_plot.jl")
include("visualization/slider_plot.jl")
include("visualization/globe.jl")
include("get_FTLE.jl")

export FTLEResult
export final_ftle
export final_ftle_field
export ftle_field
export stretching_factor
export stretching_factor!
export FTLE_from_particles!
export FTLE_from_particles
export FTLE_from_particle_file!
export FTLE_from_particle_file
export surface_plot
export slider_plot
export SliderPlotHandle
export animate_slider_plot
export globe_plot
export get_FTLE
export initial_FTLE_particle_positions!
export initial_FTLE_particle_positions
export positive_FTLE
export negative_FTLE
export Re

end
