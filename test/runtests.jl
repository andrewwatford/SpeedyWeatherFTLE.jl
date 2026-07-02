using SpeedyWeatherFTLE
using Test

include("FTLE_computations.jl")
include("get_FTLE_linear_integration.jl")
include("FTLE_result.jl")
include("get_FTLE.jl")

plotting_tests_available = try
    @eval using CairoMakie
    @eval using GeoMakie
    CairoMakie.activate!()
    true
catch err
    @info "Skipping optional plotting tests; add CairoMakie and GeoMakie to the active test environment to run them" error=sprint(showerror, err)
    false
end

if plotting_tests_available
    include("visualization/surface_plot.jl")
    include("visualization/slider_plot.jl")
    include("visualization/globe.jl")
end
