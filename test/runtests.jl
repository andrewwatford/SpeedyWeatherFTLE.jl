using NCDatasets
using RingGrids
using SpeedyWeather
using SpeedyWeatherFTLE
using Test

function stretched_particle_history(dist_km, rates, times)
    delta = rad2deg(dist_km * 1000 / SpeedyWeatherFTLE.Re)
    plonds = Matrix{Float64}(undef, 4, length(times))
    platds = similar(plonds)

    for (j, time) in enumerate(times)
        sx = exp(rates[1] * time)
        sy = exp(rates[2] * time)
        plonds[:, j] .= (sx * delta, -sx * delta, 0.0, 0.0)
        platds[:, j] .= (0.0, 0.0, sy * delta, -sy * delta)
    end

    return plonds, platds
end

@testset "SpeedyWeatherFTLE" begin
    dist_km = 10.0

    @testset "particle stencil" begin
        londs = [10.0, 20.0]
        latds = [0.0, 30.0]
        plonds, platds = initial_FTLE_particle_positions(londs, latds, dist_km)

        @test length(plonds) == 8
        @test length(platds) == 8
        @test platds[1:4] == [0.0, 0.0, platds[3], platds[4]]
        @test plonds[1] > londs[1]
        @test plonds[2] < londs[1]

        buf_lon = similar(plonds)
        buf_lat = similar(platds)
        @test initial_FTLE_particle_positions!(buf_lon, buf_lat, londs, latds, dist_km) === (buf_lon, buf_lat)
        @test buf_lon ≈ plonds
        @test buf_lat ≈ platds
        @test_throws ArgumentError initial_FTLE_particle_positions([0.0], [89.9], dist_km)
    end

    @testset "FTLE from particles" begin
        times = [0.0, 2.0, 4.0]
        rates = (0.25, -0.10)
        plonds, platds = stretched_particle_history(dist_km, rates, times)

        ftle, selected = FTLE_from_particles(plonds, platds, times, 1, dist_km; time_indices=:nonzero)

        @test selected == [2.0, 4.0]
        @test size(ftle) == (1, 2)
        @test vec(ftle) ≈ fill(rates[1], 2) rtol=1e-6

        all_ftle, all_times = FTLE_from_particles(plonds, platds, times, 1, dist_km)
        @test all_times == times
        @test isnan(all_ftle[1, 1])
        @test vec(all_ftle[:, 2:end]) ≈ vec(ftle)

        out = fill(NaN, 1, 1)
        B = fill(NaN, 2, 2, 1)
        @test FTLE_from_particles!(out, B, plonds, platds, times, 1, dist_km; time_indices=:last) === out
        @test only(out) ≈ rates[1] rtol=1e-6
        @test_throws BoundsError FTLE_from_particles(plonds, platds, times, 1, dist_km; time_indices=4)
    end

    @testset "particle files" begin
        times = [0.0, 3.0]
        rates = (0.15, 0.05)
        plonds, platds = stretched_particle_history(dist_km, rates, times)
        path = tempname() * ".nc"

        ds = NCDataset(path, "c")
        try
            defDim(ds, "particle", 4)
            defDim(ds, "time", length(times))
            defVar(ds, "time", Float64, ("time",))[:] = times
            defVar(ds, "lon", Float64, ("particle", "time"))[:, :] = plonds
            defVar(ds, "lat", Float64, ("particle", "time"))[:, :] = platds
        finally
            close(ds)
        end

        try
            ftle, selected = FTLE_from_particle_file(path, 1, dist_km; time_indices=:last, validate_initial_positions=false)
            @test selected == [3.0]
            @test only(ftle) ≈ rates[1] rtol=1e-6
        finally
            rm(path; force=true)
        end

        grid = FullClenshawGrid(4)
        plonds0, platds0 = initial_FTLE_particle_positions(grid, dist_km)
        path = tempname() * ".nc"

        ds = NCDataset(path, "c")
        try
            defDim(ds, "particle", length(plonds0))
            defDim(ds, "time", 2)
            defVar(ds, "time", Float64, ("time",))[:] = [0.0, 2.0]
            defVar(ds, "lon", Float64, ("particle", "time"))[:, :] = repeat(plonds0, 1, 2)
            defVar(ds, "lat", Float64, ("particle", "time"))[:, :] = repeat(platds0, 1, 2)
        finally
            close(ds)
        end

        try
            ftle, selected = FTLE_from_particle_file(path, grid, dist_km; time_indices=:last)
            @test selected == [2.0]
            @test all(value -> isapprox(value, 0.0; atol=1e-12), ftle)

            ds = NCDataset(path, "a")
            try
                ds.attrib["SpeedyWeatherFTLE_dist_km"] = 20.0
            finally
                close(ds)
            end
            @test_throws ArgumentError FTLE_from_particle_file(path, grid, dist_km; time_indices=:last)
        finally
            rm(path; force=true)
        end
    end

    @testset "result helpers" begin
        spectral_grid = SpectralGrid(nlayers=1, trunc=4, Grid=FullClenshawGrid)
        ftle = reshape(collect(1.0:(2 * spectral_grid.npoints)), spectral_grid.npoints, 2)
        result = FTLEResult(ftle, spectral_grid, [0.0, 2.0]; dist_km, backwards=false, rint_hours=2)

        @test size(result) == size(ftle)
        @test result[1, 2] == ftle[1, 2]
        @test final_ftle(result) == ftle[:, end]
        @test ftle_field(result; time_indices=:last) isa Field
        @test final_ftle_field(result) isa Field
        @test ftle_field(result; time_hour=1.6) isa Field
        @test stretching_factor(result)[:, 1] ≈ ones(spectral_grid.npoints)
        @test stretching_factor(ftle[:, 2], 2.0) ≈ exp.(ftle[:, 2] .* 2.0)
        @test occursin("FTLEResult", sprint(show, result))

        @test_throws DimensionMismatch FTLEResult(ftle[1:end - 1, :], spectral_grid, [0.0, 2.0]; dist_km)
        @test_throws DimensionMismatch FTLEResult(ftle, spectral_grid, [0.0]; dist_km)
        @test_throws ArgumentError FTLEResult(ftle, spectral_grid, [0.0, NaN]; dist_km)
        @test_throws ArgumentError ftle_field(FTLEResult(ftle, nothing, [0.0, 2.0]; dist_km))
    end

    @testset "frozen-field guardrails" begin
        grid = FullGaussianGrid(4)
        u = rand(grid)
        v = rand(grid)

        @test_throws ArgumentError get_FTLE(u, v; dynamics=true)
        @test_throws ArgumentError positive_FTLE(u, v; backwards=true)
        @test_throws ArgumentError negative_FTLE(u, v; backwards=false)

        result = positive_FTLE(
            u,
            v;
            simulation_days=0.25,
            rint_hours=3,
            particle_advection_every_n_time_steps=1,
            return_result=true,
            time_indices=:last,
        )

        @test result isa FTLEResult
        @test size(result, 1) == result.spectral_grid.npoints
        @test size(result, 2) == 1
        @test result.time_hours == [6.0]
        @test result.direction == :positive
        @test all(isfinite, result.ftle)
    end

    @testset "minimal dependency surface" begin
        loaded_names = Set(string(pkg.name) for pkg in keys(Base.loaded_modules))
        @test !("Makie" in loaded_names)
        @test !("GeoMakie" in loaded_names)
    end
end
