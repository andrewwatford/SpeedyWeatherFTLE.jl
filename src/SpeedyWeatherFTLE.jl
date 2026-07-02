module SpeedyWeatherFTLE

using SpeedyWeather, RingGrids
using LinearAlgebra
using NCDatasets
using Requires: @require

"""
    Re

Average Earth radius used by SpeedyWeatherFTLE, in metres.
"""
const Re = 6.371e6 # Average Earth radius in meters

include("grid_helpers.jl")
include("./FTLE_computations.jl")
include("visualization/ftle_field.jl")
include("FTLE_result.jl")
include("visualization/plot_helpers.jl")
include("plotting_api.jl")
include("FTLE_diagnostics.jl")
include("get_FTLE.jl")

function __init__()
    @require GeoMakie = "db073c08-6b98-4ee5-b6a4-5efafb3259c6" include("visualization/makie_methods.jl")
end

export FTLEResult
export final_ftle
export final_ftle_field
export ftle_field
export ftle_colorrange
export stretching_factor
export stretching_factor!
export FTLE_from_particles!
export FTLE_from_particles
export FTLE_from_particle_file!
export FTLE_from_particle_file
export surface_plot
export slider_plot
export SliderPlotHandle
export set_slider_time!
export animate_slider_plot
export globe_plot
export get_FTLE
export FTLEParticleSetup
export initial_FTLE_particle_positions!
export initial_FTLE_particle_positions
export prepare_FTLE_particles!
export attach_FTLE_tracker!
export positive_FTLE
export negative_FTLE
export Re

end
