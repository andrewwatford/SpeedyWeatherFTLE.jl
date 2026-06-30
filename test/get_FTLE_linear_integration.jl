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
end
