using SpeedyWeather, RingGrids
using SpeedyWeatherFTLE

@testset "lyapunov_FTLE.jl" begin
    # Set up a random u and v field on a FullGaussianGrid
    u_data = 1000 * randn(Float64, 8, 4)
    u_field = FullGaussianGrid(u_data, input_as=Matrix)
    v_data = 1000 * randn(Float64, 8, 4)
    v_field = FullGaussianGrid(v_data, input_as=Matrix)
    for dynamics in (true, false)
        for backwards in (true, false)
            for use_climatological in (true, false)
                use_random = !use_climatological
                gFTLE, ggrid, time_hours = lyapunov_FTLE(; dynamics=dynamics, backwards=backwards, use_climatological=use_climatological, use_random=use_random)
                @test isa(gFTLE, Matrix{Float64})
                @test isa(ggrid, AbstractGrid)
                @test isa(time_hours, Vector{Float64})
            end
            gFTLE, ggrid, time_hours = lyapunov_FTLE(; dynamics=dynamics, backwards=backwards, use_climatological=false, use_random=false, zonal_velocity_field=u_field, meridional_velocity_field=v_field)
            @test isa(gFTLE, Matrix{Float64})
            @test isa(ggrid, AbstractGrid)
            @test isa(time_hours, Vector{Float64})
        end
    end
end