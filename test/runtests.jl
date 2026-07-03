using RingGrids
using SpeedyWeather
using SpeedyWeatherFTLE
using Test

@testset "SpeedyWeatherFTLE" begin
    dist_km = 10.0

    @testset "public API surface" begin
        @test isdefined(SpeedyWeatherFTLE, :FTLE)
        @test !isdefined(SpeedyWeatherFTLE, :get_FTLE)
        @test !isdefined(SpeedyWeatherFTLE, :positive_FTLE)
        @test !isdefined(SpeedyWeatherFTLE, :negative_FTLE)
        @test !isdefined(SpeedyWeatherFTLE, :FTLE_from_particle_file)
        @test !isdefined(SpeedyWeatherFTLE, :FTLE_from_particles)
        @test !isdefined(SpeedyWeatherFTLE, :initial_FTLE_particle_positions)
    end

    @testset "result helpers" begin
        spectral_grid = SpectralGrid(nlayers=1, trunc=4, Grid=FullClenshawGrid)
        ftle = reshape(collect(1.0:(2 * spectral_grid.npoints)), spectral_grid.npoints, 2)
        result = FTLEResult(ftle, spectral_grid, [0.0, 2.0]; dist_km, rint_hours=2)
        backward = FTLEResult(ftle, spectral_grid, [0.0, 2.0]; dist_km, backwards=true, rint_hours=2)

        @test size(result) == size(ftle)
        @test result[1, 2] == ftle[1, 2]
        @test result.direction == :forward
        @test backward.direction == :backward
        @test !hasproperty(result, :particle_file_path)
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

    @testset "flow-field FTLE" begin
        grid = FullGaussianGrid(4)
        u = 0 .* rand(grid)
        v = 0 .* rand(grid)

        @test_throws ArgumentError FTLE(u, v; dynamics=true)

        result = FTLE(
            u,
            v;
            backwards=true,
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
        @test result.backwards
        @test result.direction == :backward
        @test all(isfinite, result.ftle)
        @test maximum(abs, result.ftle) < 1e-2
    end

    @testset "minimal dependency surface" begin
        loaded_names = Set(string(pkg.name) for pkg in keys(Base.loaded_modules))
        @test !("Makie" in loaded_names)
        @test !("GeoMakie" in loaded_names)
    end

    @testset "Makie plotting extension" begin
        @eval using CairoMakie
        @eval using GeoMakie
        CairoMakie.activate!()

        spectral_grid = SpectralGrid(nlayers=1, trunc=4, Grid=FullClenshawGrid)
        ftle = zeros(spectral_grid.npoints, 3)
        ftle[:, 1] .= NaN
        ftle[:, 2] .= range(0.0, 0.2; length=spectral_grid.npoints)
        ftle[:, 3] .= range(0.1, 0.3; length=spectral_grid.npoints)
        result = FTLEResult(ftle, spectral_grid, [0.0, 3.0, 6.0]; dist_km, rint_hours=3)
        lon = collect(-180:90:180)
        lat = collect(-90:45:90)

        fig, ax, sp, cb = surface_plot(result; lon, lat, coastlines=false, colorbar=false)
        @test fig !== nothing
        @test ax !== nothing
        @test sp !== nothing
        @test cb === nothing

        fig_hour, _, _, _ = surface_plot(result; time_hour=5.0, lon, lat, coastlines=false, colorbar=false)
        @test fig_hour !== nothing

        handle = slider_plot(result; lon, lat, coastlines=false, colorbar=false, return_handle=true)
        @test handle isa SliderPlotHandle
        @test handle.times == [3.0, 6.0]
        @test set_slider_time!(handle, 6.0) === handle

        function fake_record(frame_function, fig, path, frames; framerate, kwargs...)
            for frame in frames
                frame_function(frame)
            end
            return path
        end
        @test animate_slider_plot(
            "synthetic-ftle.gif",
            result;
            lon,
            lat,
            coastlines=false,
            colorbar=false,
            record_function=fake_record,
        ) == "synthetic-ftle.gif"

        globe_fig, globe_ax, globe_sp, globe_cb = globe_plot(
            result;
            lon,
            lat,
            coastlines=false,
            colorbar=false,
            show_axis=false,
        )
        @test globe_fig !== nothing
        @test globe_ax isa GeoMakie.GlobeAxis
        @test globe_sp !== nothing
        @test globe_cb === nothing
    end
end
