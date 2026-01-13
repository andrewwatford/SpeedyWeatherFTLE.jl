using SpeedyWeather, RingGrids
using SpeedyWeatherFTLE

@testset "lyapunov_FTLE.jl" begin
    for dynamics in (true, false)
        for backwards in (true, false)
            for use_climatological in (true, false)
                gFTLE, ggrid, time_hours = lyapunov_FTLE(; dynamics=dynamics, backwards=backwards, use_climatological=use_climatological)
                @test isa(gFTLE, Matrix{Float64})
                @test isa(ggrid, AbstractGrid)
                @test isa(time_hours, Vector{Float64})
            end
        end
    end
end