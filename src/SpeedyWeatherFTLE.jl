module SpeedyWeatherFTLE

using Dates: Hour, Second
using NCDatasets
using RingGrids
using SpeedyWeather

export Re
export FTLE
export shared_colorrange
export stretching_factor, stretching_factor!

export SliderPlotHandle
export animate_slider_plot
export globe_plot
export set_slider_time!
export slider_plot
export surface_plot

"""
    Re

Default planet radius used by SpeedyWeatherFTLE, in metres.

The value is inherited from SpeedyWeather.
"""
const Re = SpeedyWeather.DEFAULT_RADIUS

const _MAX_STENCIL_DEGREES = 20.0
const _PARTICLE_FILE_METADATA_PREFIX = "SpeedyWeatherFTLE_"
const _DIST_KM_ATTRIBUTE = _PARTICLE_FILE_METADATA_PREFIX * "dist_km"
const _RADIUS_ATTRIBUTE = _PARTICLE_FILE_METADATA_PREFIX * "radius"
const _PARTICLE_ORDER_ATTRIBUTE = _PARTICLE_FILE_METADATA_PREFIX * "particle_order"
const _PARTICLE_ORDER = "east,west,north,south"

include("utilities.jl")
include("particles.jl")
include("FTLE_computations.jl")
include("FTLE_diagnostics.jl")
include("FTLE.jl")
include("plotting_api.jl")

end
