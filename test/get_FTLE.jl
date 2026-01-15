using SpeedyWeather, RingGrids
using SpeedyWeatherFTLE
using InteractiveUtils

@testset "get_FTLE.jl" begin
    nlat_half = 5
    full_grids = subtypes(RingGrids.AbstractFullGrid)
    reduced_grids = subtypes(RingGrids.AbstractReducedGrid)
    all_grids = vcat(full_grids, reduced_grids)
    for grid in all_grids
        spatial_grid = grid(nlat_half)
        u = rand(spatial_grid)
        v = rand(spatial_grid)

        for dynamics in (true, false)
            for backwards in (true, false)
                FTLE, spectral_grid, time_hours = get_FTLE(u, v; dynamics=dynamics, backwards=backwards)
                @test isa(FTLE, Matrix{Float64})
                @test isa(spectral_grid, SpectralGrid)
                @test isa(time_hours, Vector{Float64})
            end
        end
    end
end
