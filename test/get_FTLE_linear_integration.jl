using LinearAlgebra
using Logging
using RingGrids
using SpeedyWeather
using SpeedyWeatherFTLE
using Test

@testset "get_FTLE linear integration" begin
    center_lond = 180.0
    grid = FullClenshawGrid(4)

    wrapped_lond_delta(lond, center_lond) = mod(lond - center_lond + 180, 360) - 180

    function linear_velocity_fields(A)
        # Build a local tangent-plane linear velocity field u = A*x around the
        # equatorial center point. The cos(latitude) factor converts eastward
        # physical velocity into the spherical zonal velocity convention used by
        # SpeedyWeather's particle advection.
        londs, latds = RingGrids.get_londlatds(grid)
        x = SpeedyWeatherFTLE.Re .* deg2rad.(wrapped_lond_delta.(londs, center_lond))
        y = SpeedyWeatherFTLE.Re .* deg2rad.(latds)

        u = rand(grid)
        v = rand(grid)
        u .= cosd.(latds) .* (A[1, 1] .* x .+ A[1, 2] .* y)
        v .= A[2, 1] .* x .+ A[2, 2] .* y

        return u, v
    end

    function center_gridpoint(grid)
        londs, latds = RingGrids.get_londlatds(grid)
        return only(findall(i -> londs[i] == center_lond && latds[i] == 0, eachindex(londs)))
    end

    # These cases cover stretching, contraction, shear, non-normal dynamics, and
    # pure rotation. Tolerances include SpeedyWeather interpolation and Heun-step
    # error; the center-point cases remain tight for diagonal and simple shear flows.
    linear_systems = [
        ("zonal strain", [2e-6 0.0; 0.0 0.0], 1e-2, 2e-5),
        ("zonal contraction", [-2e-6 0.0; 0.0 0.0], 1e-2, 2e-5),
        ("meridional strain", [0.0 0.0; 0.0 2e-6], 1e-2, 2e-5),
        ("meridional shear", [0.0 0.0; 1e-6 0.0], 1e-2, 2e-5),
        ("mixed shear", [0.0 1e-6; -5e-7 0.0], 2e-1, 2e-4),
        ("pure rotation", [0.0 -1e-6; 1e-6 0.0], 1e-2, 2e-4),
    ]

    for (_, A, rtol, atol) in linear_systems
        for backwards in (false, true)
            u, v = linear_velocity_fields(A)
            FTLE, spectral_grid, time_hours = with_logger(ConsoleLogger(stderr, Logging.Warn)) do
                get_FTLE(
                    u,
                    v;
                    simulation_days = 0.75,
                    dist_km = 500,
                    dynamics = false,
                    backwards,
                    rint_hours = 6,
                    particle_tracker_keepbits = 23,
                )
            end

            # FullClenshawGrid(4) advances particles every 18 hours, so ending
            # the simulation at 0.75 days compares against an actual advection update.
            t = time_hours[end]
            @test t == 18

            time_direction = backwards ? -1 : 1
            flow_map = exp(time_direction * 3600 * t * A)
            expected = log(maximum(svdvals(flow_map))) / t
            actual = FTLE[center_gridpoint(spectral_grid.grid), end]

            @test actual ≈ expected rtol=rtol atol=atol
        end
    end

    @testset "particle advection cadence keyword" begin
        # This runs the same linear-flow check with particles updated every model
        # time step. The recorded final time should now be a genuine 6-hour
        # advection result, rather than waiting for the default 18-hour cadence.
        A = [2e-6 0.0; 0.0 0.0]
        u, v = linear_velocity_fields(A)
        tracker_dir = mktempdir()
        particle_path = nothing

        try
            FTLE, spectral_grid, time_hours, particle_path = with_logger(ConsoleLogger(stderr, Logging.Warn)) do
                get_FTLE(
                    u,
                    v;
                    simulation_days = 0.25,
                    dist_km = 500,
                    dynamics = false,
                    backwards = false,
                    rint_hours = 3,
                    particle_advection_every_n_time_steps = 1,
                    particle_tracker_keepbits = 23,
                    particle_tracker_path = tracker_dir,
                    particle_tracker_filename = "linear_particles.nc",
                    return_particle_file_path = true,
                )
            end

            @test isfile(particle_path)

            reread_FTLE, reread_time_hours = FTLE_from_particle_file(particle_path, spectral_grid, 500)
            @test isequal(reread_FTLE, FTLE)
            @test reread_time_hours == time_hours

            FTLE_buffer = fill(NaN, size(FTLE))
            B_buffer = fill(NaN, 2, 2, spectral_grid.npoints)
            inplace_FTLE, inplace_time_hours = FTLE_from_particle_file!(
                FTLE_buffer,
                B_buffer,
                particle_path,
                spectral_grid,
                500,
            )

            @test inplace_FTLE === FTLE_buffer
            @test isequal(inplace_FTLE, FTLE)
            @test inplace_time_hours == time_hours

            subset_indices = [length(time_hours), length(time_hours) - 1]
            subset_FTLE, subset_time_hours = FTLE_from_particle_file(
                particle_path,
                spectral_grid,
                500;
                time_indices = subset_indices,
            )
            @test isequal(subset_FTLE, FTLE[:, subset_indices])
            @test subset_time_hours == time_hours[subset_indices]

            nonzero_FTLE, nonzero_time_hours = FTLE_from_particle_file(
                particle_path,
                spectral_grid,
                500;
                time_indices = :nonzero,
            )
            @test isequal(nonzero_FTLE, FTLE[:, 2:end])
            @test nonzero_time_hours == time_hours[2:end]

            single_time_FTLE = fill(NaN, spectral_grid.npoints, 1)
            single_time_result, single_time_hours = FTLE_from_particle_file!(
                single_time_FTLE,
                B_buffer,
                particle_path,
                spectral_grid,
                500;
                time_indices = :last,
            )
            @test single_time_result === single_time_FTLE
            @test isequal(single_time_result, FTLE[:, end:end])
            @test single_time_hours == [time_hours[end]]

            @test_throws DimensionMismatch FTLE_from_particle_file!(
                fill(NaN, spectral_grid.npoints + 1, length(time_hours)),
                B_buffer,
                particle_path,
                spectral_grid,
                500,
            )
            @test_throws DimensionMismatch FTLE_from_particle_file!(
                FTLE_buffer,
                fill(NaN, 2, 2, spectral_grid.npoints + 1),
                particle_path,
                spectral_grid,
                500,
            )
            @test_throws BoundsError FTLE_from_particle_file(particle_path, spectral_grid, 500; time_indices = [length(time_hours) + 1])
            @test_throws ArgumentError FTLE_from_particle_file(particle_path, spectral_grid, 500; time_indices = 1.5)

            t = time_hours[end]
            @test t == 6

            flow_map = exp(3600 * t * A)
            expected = log(maximum(svdvals(flow_map))) / t
            actual = FTLE[center_gridpoint(spectral_grid.grid), end]

            @test actual ≈ expected rtol=1e-2 atol=2e-5
        finally
            particle_path === nothing || rm(particle_path; force=true)
            rm(tracker_dir; force=true, recursive=true)
        end

        @test_throws ArgumentError get_FTLE(u, v; particle_advection_every_n_time_steps = 0)
        @test_throws ArgumentError get_FTLE(u, v; dist_km = 0)
        @test_throws ArgumentError get_FTLE(u, v; particle_tracker_compression_level = 10)
    end
end
