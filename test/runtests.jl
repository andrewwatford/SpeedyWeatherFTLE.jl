using RingGrids
using SpeedyWeather
using SpeedyWeatherFTLE
using Test

@testset "SpeedyWeatherFTLE" begin
    @testset "public API surface" begin
        @test isdefined(SpeedyWeatherFTLE, :FTLE)
        @test isdefined(SpeedyWeatherFTLE, :shared_colorrange)
        @test !isdefined(SpeedyWeatherFTLE, :FTLEResult)
        @test !isdefined(SpeedyWeatherFTLE, :final_ftle)
        @test !isdefined(SpeedyWeatherFTLE, :final_ftle_field)
        @test !isdefined(SpeedyWeatherFTLE, :ftle_field)
        @test !isdefined(SpeedyWeatherFTLE, :ftle_colorrange)
        @test !isdefined(SpeedyWeatherFTLE, :get_FTLE)
        @test !isdefined(SpeedyWeatherFTLE, :positive_FTLE)
        @test !isdefined(SpeedyWeatherFTLE, :negative_FTLE)
        @test !isdefined(SpeedyWeatherFTLE, :FTLE_from_particle_file)
        @test !isdefined(SpeedyWeatherFTLE, :FTLE_from_particles)
        @test !isdefined(SpeedyWeatherFTLE, :initial_FTLE_particle_positions)
    end

    @testset "field diagnostics" begin
        spectral_grid = SpectralGrid(nlayers=1, trunc=4, Grid=FullClenshawGrid)
        values = reshape(collect(1.0:(2 * spectral_grid.npoints)), spectral_grid.npoints, 2)
        field = Field(values, spectral_grid.grid)

        @test shared_colorrange(field) == (1.0, maximum(values))
        @test shared_colorrange([-2.0, 1.0]; symmetric=true) == (-2.0, 2.0)

        stretch = stretching_factor(field, [0.0, 2.0])
        @test stretch isa Field
        @test stretch[:, 1] ≈ Field(ones(spectral_grid.npoints), spectral_grid.grid)
        @test stretch[:, 2] ≈ Field(exp.(values[:, 2] .* 2.0), spectral_grid.grid)
        @test stretching_factor(field[:, 2], 2.0) ≈ stretch[:, 2]
    end

    @testset "flow-field FTLE" begin
        grid = FullGaussianGrid(4)
        u = 0 .* rand(grid)
        v = 0 .* rand(grid)

        @test_throws MethodError FTLE(u, v; dynamics=true)

        ftle, time_hours = FTLE(
            u,
            v;
            backwards=true,
            simulation_days=0.25,
            rint_hours=3,
            particle_advection_every_n_time_steps=1,
            time_indices=:last,
        )

        @test ftle isa Field
        @test size(ftle, 1) == length(first(RingGrids.get_londlatds(ftle.grid)))
        @test size(ftle, 2) == 1
        @test time_hours == [6.0]
        @test all(isfinite, ftle)
        @test maximum(abs, ftle) < 1e-2
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
        values = zeros(spectral_grid.npoints, 3)
        values[:, 1] .= NaN
        values[:, 2] .= range(0.0, 0.2; length=spectral_grid.npoints)
        values[:, 3] .= range(0.1, 0.3; length=spectral_grid.npoints)
        ftle = Field(values, spectral_grid.grid)
        time_hours = [0.0, 3.0, 6.0]
        lon = collect(-180:90:180)
        lat = collect(-90:45:90)

        fig, ax, sp, cb = surface_plot(ftle; time_hours, lon, lat, coastlines=false, colorbar=false)
        @test fig !== nothing
        @test ax !== nothing
        @test sp !== nothing
        @test cb === nothing

        fig_hour, _, _, _ = surface_plot(ftle; time_hours, time_hour=5.0, lon, lat, coastlines=false, colorbar=false)
        @test fig_hour !== nothing

        handle = slider_plot(time_hours, ftle; lon, lat, coastlines=false, colorbar=false, return_handle=true)
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
            time_hours,
            ftle;
            lon,
            lat,
            coastlines=false,
            colorbar=false,
            record_function=fake_record,
        ) == "synthetic-ftle.gif"

        globe_fig, globe_ax, globe_sp, globe_cb = globe_plot(
            ftle;
            time_hours,
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
