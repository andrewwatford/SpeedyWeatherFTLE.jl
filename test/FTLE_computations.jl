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
        actual = fill(NaN, length(linear_systems))

        @test SpeedyWeatherFTLE.FTLE_over_grid!(actual, B, duration) === actual
        @test actual ≈ expected rtol=1e-12 atol=1e-12
        @test SpeedyWeatherFTLE.FTLE_over_grid(B, duration) ≈ expected rtol=1e-12 atol=1e-12
    end

    @testset "negative-time linear systems" begin
        # Negative-time FTLE uses the backward flow map exp(-A*T), which can differ
        # from the positive-time value for contracting or non-normal systems.
        B = exact_flow_maps(linear_systems, duration, -1)
        expected = [exact_ftle(A, duration, -1) for (_, A) in linear_systems]

        @test SpeedyWeatherFTLE.FTLE_over_grid(B, duration) ≈ expected rtol=1e-12 atol=1e-12
    end

    @testset "zero-duration samples are undefined" begin
        B = exact_flow_maps(linear_systems, 0.0, 1)
        actual = fill(0.0, length(linear_systems))

        @test SpeedyWeatherFTLE.FTLE_over_grid!(actual, B, 0.0) === actual
        @test all(isnan, actual)
        @test all(isnan, SpeedyWeatherFTLE.FTLE_over_grid(B, 0.0))
    end

    @testset "wrapped longitude displacement gradient" begin
        # East-west particles can straddle 0/360 degrees; their central difference
        # should still reconstruct an identity deformation at release time.
        dist_km = 10.0
        delta_degrees = rad2deg(dist_km * 1000 / SpeedyWeatherFTLE.Re)
        plonds = [delta_degrees, 360 - delta_degrees, 0.0, 0.0]
        platds = [0.0, 0.0, delta_degrees, -delta_degrees]
        B = fill(NaN, 2, 2, 1)

        @test SpeedyWeatherFTLE.displacement_gradient_matrix_central!(B, plonds, platds, dist_km) === B
        @test B[:, :, 1] ≈ I
        @test SpeedyWeatherFTLE.displacement_gradient_matrix_central(plonds, platds, dist_km)[:, :, 1] ≈ I
    end

    @testset "initial FTLE particle release positions" begin
        dist_km = 10.0
        londs = [10.0, 20.0]
        latds = [0.0, 30.0]
        delta_degrees = rad2deg(dist_km * 1000 / SpeedyWeatherFTLE.Re)

        plonds, platds = initial_FTLE_particle_positions(londs, latds, dist_km)
        expected_plonds = [
            londs[1] + delta_degrees,
            londs[1] - delta_degrees,
            londs[1],
            londs[1],
            londs[2] + delta_degrees / cosd(latds[2]),
            londs[2] - delta_degrees / cosd(latds[2]),
            londs[2],
            londs[2],
        ]
        expected_platds = [
            latds[1],
            latds[1],
            latds[1] + delta_degrees,
            latds[1] - delta_degrees,
            latds[2],
            latds[2],
            latds[2] + delta_degrees,
            latds[2] - delta_degrees,
        ]

        @test plonds ≈ expected_plonds
        @test platds ≈ expected_platds

        plonds_buffer = fill(NaN, length(plonds))
        platds_buffer = fill(NaN, length(platds))
        returned_plonds, returned_platds = initial_FTLE_particle_positions!(
            plonds_buffer,
            platds_buffer,
            londs,
            latds,
            dist_km,
        )
        @test returned_plonds === plonds_buffer
        @test returned_platds === platds_buffer
        @test plonds_buffer ≈ plonds
        @test platds_buffer ≈ platds

        particles = fill(Particle(0.0, 0.0), length(plonds))
        SpeedyWeatherFTLE.perturb_positions_FTLE(particles, londs, latds, dist_km)
        @test [particle.lon for particle in particles] ≈ plonds
        @test [particle.lat for particle in particles] ≈ platds

        @test_throws DimensionMismatch initial_FTLE_particle_positions(londs, latds[1:1], dist_km)
        @test_throws DimensionMismatch initial_FTLE_particle_positions!(plonds_buffer[1:end - 1], platds_buffer, londs, latds, dist_km)
        @test_throws DimensionMismatch initial_FTLE_particle_positions!(plonds_buffer, platds_buffer[1:end - 1], londs, latds, dist_km)
        @test_throws DimensionMismatch SpeedyWeatherFTLE.perturb_positions_FTLE(particles[1:end - 1], londs, latds, dist_km)
        @test_throws ArgumentError initial_FTLE_particle_positions(londs, latds, 0)
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

        function linear_trajectory_positions(A, time_direction, times)
            particles = fill(Particle(center_lond, 0.0), 4)
            SpeedyWeatherFTLE.perturb_positions_FTLE(particles, londs, latds, dist_km)

            plonds_time = Matrix{Float64}(undef, length(particles), length(times))
            platds_time = similar(plonds_time)

            for (tindex, time) in enumerate(times)
                flow_map = exp(time_direction * time * A)
                for (pindex, particle) in enumerate(particles)
                    initial_position = SpeedyWeatherFTLE.Re .* [
                        deg2rad(particle.lon - center_lond),
                        deg2rad(particle.lat),
                    ]
                    final_position = flow_map * initial_position
                    plonds_time[pindex, tindex] = center_lond + rad2deg(final_position[1] / SpeedyWeatherFTLE.Re)
                    platds_time[pindex, tindex] = rad2deg(final_position[2] / SpeedyWeatherFTLE.Re)
                end
            end

            return plonds_time, platds_time
        end

        for (_, A) in linear_systems
            expected_positive = exact_ftle(A, duration, 1)
            expected_negative = exact_ftle(A, duration, -1)

            @test linear_trajectory_ftle(A, 1) ≈ expected_positive rtol=1e-4 atol=1e-8
            @test linear_trajectory_ftle(A, -1) ≈ expected_negative rtol=1e-4 atol=1e-8
        end

        @testset "FTLE_from_particles" begin
            times = [0.0, duration, 2duration]
            A = [0.05 1.20; -0.40 -0.15]
            plonds_time, platds_time = linear_trajectory_positions(A, 1, times)
            FTLE_grid_time, returned_times = FTLE_from_particles(plonds_time, platds_time, times, 1, dist_km)
            FTLE_buffer = fill(NaN, 1, length(times))
            B_buffer = fill(NaN, 2, 2, 1)
            expected = [exact_ftle(A, time, 1) for time in times[2:end]]

            @test returned_times == Float64.(times)
            @test all(isnan, FTLE_grid_time[:, 1])
            @test vec(FTLE_grid_time[:, 2:end]) ≈ expected rtol=1e-4 atol=1e-8

            expected_stretching = exp.(expected .* times[2:end])
            stretching = stretching_factor(FTLE_grid_time, returned_times)
            stretching_buffer = fill(NaN, size(FTLE_grid_time))
            @test stretching_factor!(stretching_buffer, FTLE_grid_time, returned_times) === stretching_buffer
            @test isequal(stretching_buffer, stretching)
            @test all(isnan, stretching[:, 1])
            @test vec(stretching[:, 2:end]) ≈ expected_stretching rtol=1e-4 atol=1e-8
            @test stretching_factor(vec(FTLE_grid_time[:, 2]), returned_times[2]) ≈ [expected_stretching[1]] rtol=1e-4 atol=1e-8
            @test isequal(
                stretching_factor(FTLEResult(
                    FTLE_grid_time,
                    nothing,
                    returned_times;
                    dist_km,
                    backwards = false,
                    dynamics = false,
                    rint_hours = 1,
                )),
                stretching,
            )
            @test_throws DimensionMismatch stretching_factor(FTLE_grid_time, returned_times[1:2])
            @test_throws DimensionMismatch stretching_factor!(fill(NaN, 1, 1), FTLE_grid_time, returned_times)

            @test FTLE_from_particles!(FTLE_buffer, B_buffer, plonds_time, platds_time, times, 1, dist_km) === FTLE_buffer
            @test all(isnan, FTLE_buffer[:, 1])
            @test FTLE_buffer[:, 2:end] ≈ FTLE_grid_time[:, 2:end] rtol=1e-4 atol=1e-8

            subset_FTLE, subset_times = FTLE_from_particles(
                plonds_time,
                platds_time,
                times,
                1,
                dist_km;
                time_indices = [3, 2],
            )
            @test subset_times == Float64.(times[[3, 2]])
            @test subset_FTLE[:, 1] ≈ FTLE_grid_time[:, 3] rtol=1e-4 atol=1e-8
            @test subset_FTLE[:, 2] ≈ FTLE_grid_time[:, 2] rtol=1e-4 atol=1e-8

            last_FTLE, last_times = FTLE_from_particles(
                plonds_time,
                platds_time,
                times,
                1,
                dist_km;
                time_indices = :last,
            )
            @test last_times == [Float64(times[end])]
            @test last_FTLE ≈ FTLE_grid_time[:, end:end] rtol=1e-4 atol=1e-8

            nonzero_FTLE, nonzero_times = FTLE_from_particles(
                plonds_time,
                platds_time,
                times,
                1,
                dist_km;
                time_indices = :nonzero,
            )
            @test nonzero_times == Float64.(times[2:end])
            @test nonzero_FTLE ≈ FTLE_grid_time[:, 2:end] rtol=1e-4 atol=1e-8

            mask_FTLE, mask_times = FTLE_from_particles(
                plonds_time,
                platds_time,
                times,
                1,
                dist_km;
                time_indices = [false, true, false],
            )
            @test mask_times == [Float64(times[2])]
            @test mask_FTLE ≈ FTLE_grid_time[:, 2:2] rtol=1e-4 atol=1e-8

            single_buffer = fill(NaN, 1, 1)
            @test FTLE_from_particles!(
                single_buffer,
                B_buffer,
                plonds_time,
                platds_time,
                times,
                1,
                dist_km;
                time_indices = 2,
            ) === single_buffer
            @test vec(single_buffer) ≈ FTLE_grid_time[:, 2] rtol=1e-4 atol=1e-8

            @test_throws DimensionMismatch FTLE_from_particles(plonds_time[1:3, :], platds_time[1:3, :], times, 1, dist_km)
            @test_throws DimensionMismatch FTLE_from_particles(plonds_time, platds_time[:, 1:1], times, 1, dist_km)
            @test_throws DimensionMismatch FTLE_from_particles(plonds_time, platds_time, times[1:1], 1, dist_km)
            @test_throws ArgumentError FTLE_from_particles(plonds_time, platds_time, times, 1, 0)
            @test_throws DimensionMismatch FTLE_from_particles!(fill(NaN, 2, length(times)), B_buffer, plonds_time, platds_time, times, 1, dist_km)
            @test_throws DimensionMismatch FTLE_from_particles!(FTLE_buffer, fill(NaN, 2, 2, 2), plonds_time, platds_time, times, 1, dist_km)
            @test_throws BoundsError FTLE_from_particles(plonds_time, platds_time, times, 1, dist_km; time_indices = [4])
            @test_throws ArgumentError FTLE_from_particles(plonds_time, platds_time, times, 1, dist_km; time_indices = 2.5)
            @test_throws ArgumentError FTLE_from_particles(plonds_time, platds_time, times, 1, dist_km; time_indices = :middle)
            @test_throws DimensionMismatch FTLE_from_particles(plonds_time, platds_time, times, 1, dist_km; time_indices = [true, false])
        end
    end
end
