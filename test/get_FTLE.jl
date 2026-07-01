using SpeedyWeather, RingGrids
using SpeedyWeatherFTLE
using InteractiveUtils

@testset "get_FTLE.jl" begin
    nlat_half = 10
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

    @testset "FTLEResult return" begin
        spatial_grid = FullClenshawGrid(4)
        u = rand(spatial_grid)
        v = rand(spatial_grid)

        function capture_stderr(f)
            path = tempname()
            try
                result = open(path, "w") do io
                    redirect_stderr(io) do
                        f()
                    end
                end
                return result, read(path, String)
            finally
                rm(path; force=true)
            end
        end

        result, positive_stderr = capture_stderr() do
            positive_FTLE(
                u,
                v;
                simulation_days = 0.25,
                dynamics = false,
                rint_hours = 3,
                particle_advection_every_n_time_steps = 1,
                return_result = true,
                time_indices = :nonzero,
            )
        end

        @test isa(result, FTLEResult)
        @test !occursin("Time step changed", positive_stderr)
        @test result.particle_file_path === nothing
        @test result.dist_km == 10
        @test result.backwards == false
        @test result.direction == :positive
        @test result.dynamics == false
        @test result.rint_hours == 3
        @test result.time_hours == [3.0, 6.0]
        @test size(result, 1) == result.spectral_grid.npoints
        @test size(result, 2) == 2
        @test isa(ftle_field(result), Field)
        @test length(final_ftle(result)) == result.spectral_grid.npoints
        @test isa(final_ftle_field(result), Field)
        @test occursin("FTLEResult(", sprint(show, result))

        result_display = sprint(show, MIME"text/plain"(), result)
        @test occursin("grid points: $(result.spectral_grid.npoints)", result_display)
        @test occursin("times: 2 (3.0 to 6.0 h)", result_display)
        @test occursin("direction: positive", result_display)
        @test !occursin("ftle:", result_display)

        negative_result, negative_stderr = capture_stderr() do
            negative_FTLE(
                u,
                v;
                simulation_days = 0.25,
                dynamics = false,
                rint_hours = 3,
                particle_advection_every_n_time_steps = 1,
                return_result = true,
                time_indices = :last,
            )
        end

        @test isa(negative_result, FTLEResult)
        @test !occursin("Time step changed", negative_stderr)
        @test negative_result.backwards == true
        @test negative_result.direction == :negative
        @test negative_result.time_hours == [6.0]
        @test size(negative_result, 1) == negative_result.spectral_grid.npoints
        @test size(negative_result, 2) == 1
        @test_throws ArgumentError positive_FTLE(u, v; backwards = true)
        @test_throws ArgumentError negative_FTLE(u, v; backwards = false)
        selected_FTLE, selected_spectral_grid, selected_time_hours = get_FTLE(
            u,
            v;
            simulation_days = 0.25,
            dynamics = false,
            rint_hours = 3,
            particle_advection_every_n_time_steps = 1,
            time_indices = [false, true, true],
        )

        @test isa(selected_FTLE, Matrix{Float64})
        @test isa(selected_spectral_grid, SpectralGrid)
        @test selected_time_hours == [3.0, 6.0]
        @test size(selected_FTLE) == (selected_spectral_grid.npoints, 2)
        @test_throws ArgumentError positive_FTLE(
            u,
            v;
            simulation_days = 0.25,
            dynamics = false,
            rint_hours = 3,
            particle_advection_every_n_time_steps = 1,
            time_indices = :middle,
        )
    end
end
