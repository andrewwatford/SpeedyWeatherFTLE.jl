using LinearAlgebra
using SpeedyWeather
using SpeedyWeatherFTLE
using Test

@testset "FTLE_computations.jl" begin
    duration = 2.75

    linear_systems = [
        ("axis-aligned strain", [0.20 0.00; 0.00 -0.10]),
        ("axis-aligned contraction", [-0.15 0.00; 0.00 -0.25]),
        ("simple shear", [0.00 0.45; 0.00 0.00]),
        ("rotation with expansion", [0.08 -0.70; 0.70 0.08]),
        ("non-normal mixed flow", [0.05 1.20; -0.40 -0.15]),
    ]

    function exact_ftle(A, duration, time_direction)
        flow_map = exp(time_direction * duration * A)
        largest_singular_value = maximum(svdvals(flow_map))
        return log(largest_singular_value) / duration
    end

    function exact_flow_maps(linear_systems, duration, time_direction)
        B = Array{Float64}(undef, 2, 2, length(linear_systems))
        for (i, (_, A)) in enumerate(linear_systems)
            B[:, :, i] .= exp(time_direction * duration * A)
        end
        return B
    end

    @testset "positive-time linear systems" begin
        # If the displacement gradient is the exact forward flow map exp(A*T),
        # FTLE_over_grid should return log of the largest singular value per unit time.
        B = exact_flow_maps(linear_systems, duration, 1)
        expected = [exact_ftle(A, duration, 1) for (_, A) in linear_systems]

        @test SpeedyWeatherFTLE.FTLE_over_grid(B, duration) ≈ expected rtol=1e-12 atol=1e-12
    end

    @testset "negative-time linear systems" begin
        # Negative-time FTLE uses the backward flow map exp(-A*T), which can differ
        # from the positive-time value for contracting or non-normal systems.
        B = exact_flow_maps(linear_systems, duration, -1)
        expected = [exact_ftle(A, duration, -1) for (_, A) in linear_systems]

        @test SpeedyWeatherFTLE.FTLE_over_grid(B, duration) ≈ expected rtol=1e-12 atol=1e-12
    end

    @testset "wrapped longitude displacement gradient" begin
        # East-west particles can straddle 0/360 degrees; their central difference
        # should still reconstruct an identity deformation at release time.
        dist_km = 10.0
        delta_degrees = rad2deg(dist_km * 1000 / SpeedyWeatherFTLE.Re)
        plonds = [delta_degrees, 360 - delta_degrees, 0.0, 0.0]
        platds = [0.0, 0.0, delta_degrees, -delta_degrees]

        @test SpeedyWeatherFTLE.displacement_gradient_matrix_central(plonds, platds, dist_km)[:, :, 1] ≈ I
    end

    @testset "linear particle trajectories" begin
        # This exercises the FTLE particle layout and displacement-gradient
        # reconstruction using exact trajectories from u = A*x, without the
        # additional interpolation and time-stepping error from SpeedyWeather.
        center_lond = 180.0
        dist_km = 10.0
        londs = [center_lond]
        latds = [0.0]

        function linear_trajectory_ftle(A, time_direction)
            particles = fill(Particle(center_lond, 0.0), 4)
            SpeedyWeatherFTLE.perturb_positions_FTLE(particles, londs, latds, dist_km)

            flow_map = exp(time_direction * duration * A)
            plonds = Vector{Float64}(undef, length(particles))
            platds = Vector{Float64}(undef, length(particles))

            for (i, particle) in enumerate(particles)
                initial_position = SpeedyWeatherFTLE.Re .* [
                    deg2rad(particle.lon - center_lond),
                    deg2rad(particle.lat),
                ]
                final_position = flow_map * initial_position
                plonds[i] = center_lond + rad2deg(final_position[1] / SpeedyWeatherFTLE.Re)
                platds[i] = rad2deg(final_position[2] / SpeedyWeatherFTLE.Re)
            end

            B = SpeedyWeatherFTLE.displacement_gradient_matrix_central(plonds, platds, dist_km)
            return only(SpeedyWeatherFTLE.FTLE_over_grid(B, duration))
        end

        for (_, A) in linear_systems
            expected_positive = exact_ftle(A, duration, 1)
            expected_negative = exact_ftle(A, duration, -1)

            @test linear_trajectory_ftle(A, 1) ≈ expected_positive rtol=1e-4 atol=1e-8
            @test linear_trajectory_ftle(A, -1) ≈ expected_negative rtol=1e-4 atol=1e-8
        end
    end
end
