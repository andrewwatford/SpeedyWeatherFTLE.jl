using GeoMakie
using SpeedyWeather
using SpeedyWeatherFTLE

@testset "globe_plot.jl" begin
    grid = HEALPixGrid(20)
    field = rand(grid)
    title = "Test Globe Plot"
    label = "Field Value"

    for coastlines in (true, false)
        for colorbar in (true, false)
            fig, ax, sp, cb = globe_plot(
                field;
                lon = collect(-180:30:180),
                lat = collect(-90:30:90),
                title,
                coastlines,
                colorbar,
                label,
            )

            @test isa(fig, Figure)
            @test isa(ax, GeoMakie.GlobeAxis)
            @test isa(sp, GeoMakie.Makie.Surface)
            if colorbar
                @test isa(cb, Colorbar)
                @test cb.label[] == label
            else
                @test cb === nothing
            end
            @test ax.title[] == title
        end
    end

    fig, ax, sp, cb = globe_plot(
        field;
        lon = collect(-180:45:180),
        lat = collect(-90:45:90),
        coastlines = false,
        colorbar = true,
        figure_kwargs = (; size = (480, 320)),
        axis_kwargs = (; show_axis = true),
        surface_kwargs = (; transparency = false),
        colorbar_kwargs = (; vertical = true),
    )

    @test isa(fig, Figure)
    @test isa(ax, GeoMakie.GlobeAxis)
    @test isa(sp, GeoMakie.Makie.Surface)
    @test isa(cb, Colorbar)
    @test Tuple(fig.scene.viewport[].widths) == (480, 320)
    @test ax.show_axis[] == true

    @testset "FTLE overloads" begin
        spectral_grid = SpectralGrid(nlayers=1, trunc=6, Grid=FullGaussianGrid)
        FTLE = rand(spectral_grid.npoints, 4)
        result = FTLEResult(
            FTLE,
            spectral_grid,
            collect(0.0:3.0);
            dist_km = 10,
            backwards = false,
            dynamics = false,
            rint_hours = 1,
        )
        diagnostic_result = FTLEResult(
            FTLE,
            nothing,
            collect(0.0:3.0);
            dist_km = 10,
            backwards = false,
            dynamics = false,
            rint_hours = 1,
        )

        fig, ax, sp, cb = globe_plot(
            FTLE,
            spectral_grid;
            time_index = 2,
            lon = collect(-180:45:180),
            lat = collect(-90:45:90),
            colorbar = false,
            coastlines = false,
        )

        @test isa(fig, Figure)
        @test isa(ax, GeoMakie.GlobeAxis)
        @test isa(sp, GeoMakie.Makie.Surface)
        @test cb === nothing
        @test_throws BoundsError globe_plot(FTLE, spectral_grid; time_index = 0)
        @test_throws ArgumentError globe_plot(FTLE, spectral_grid; time_index = 2, time_hour = 2.0, time_hours = result.time_hours)
        @test_throws ArgumentError globe_plot(FTLE, spectral_grid; time_hour = 2.0)
        @test_throws ArgumentError globe_plot(FTLE, spectral_grid; time_hour = -Inf, time_hours = result.time_hours)
        @test_throws DimensionMismatch globe_plot(FTLE, spectral_grid; time_hour = 2.0, time_hours = [0.0, 1.0])

        fig, ax, sp, cb = globe_plot(
            FTLE,
            spectral_grid;
            time_hours = result.time_hours,
            time_hour = 2.6,
            lon = collect(-180:45:180),
            lat = collect(-90:45:90),
            colorbar = true,
            coastlines = false,
        )

        @test isa(fig, Figure)
        @test isa(ax, GeoMakie.GlobeAxis)
        @test isa(sp, GeoMakie.Makie.Surface)
        @test isa(cb, Colorbar)
        @test sp.colorrange[] ≈ collect(ftle_colorrange(view(FTLE, :, 4)))

        fig, ax, sp, cb = globe_plot(
            FTLE[:, end],
            spectral_grid;
            lon = collect(-180:45:180),
            lat = collect(-90:45:90),
            colorbar = true,
            coastlines = false,
        )

        @test isa(fig, Figure)
        @test isa(ax, GeoMakie.GlobeAxis)
        @test isa(sp, GeoMakie.Makie.Surface)
        @test isa(cb, Colorbar)
        @test cb.label[] == "FTLE [1/h]"
        @test sp.colorrange[] ≈ collect(ftle_colorrange(FTLE[:, end]))

        fig, ax, sp, cb = globe_plot(
            FTLE[:, end],
            spectral_grid;
            lon = collect(-180:45:180),
            lat = collect(-90:45:90),
            colorbar = false,
            coastlines = false,
        )

        @test isa(fig, Figure)
        @test isa(ax, GeoMakie.GlobeAxis)
        @test isa(sp, GeoMakie.Makie.Surface)
        @test cb === nothing
        @test_throws ArgumentError globe_plot(diagnostic_result)

        fig, ax, sp, cb = globe_plot(
            result;
            time_index = 3,
            lon = collect(-180:45:180),
            lat = collect(-90:45:90),
            colorbar = false,
            coastlines = false,
        )

        @test isa(fig, Figure)
        @test isa(ax, GeoMakie.GlobeAxis)
        @test isa(sp, GeoMakie.Makie.Surface)
        @test cb === nothing

        fig, ax, sp, cb = globe_plot(
            result;
            time_hour = 2.6,
            lon = collect(-180:45:180),
            lat = collect(-90:45:90),
            colorbar = true,
            coastlines = false,
        )

        @test isa(fig, Figure)
        @test isa(ax, GeoMakie.GlobeAxis)
        @test isa(sp, GeoMakie.Makie.Surface)
        @test isa(cb, Colorbar)
        @test sp.colorrange[] ≈ collect(ftle_colorrange(view(FTLE, :, 4)))
        @test_throws ArgumentError globe_plot(result; time_index = 2, time_hour = 2.0)
        @test_throws ArgumentError globe_plot(result; time_hour = NaN)
    end
end
