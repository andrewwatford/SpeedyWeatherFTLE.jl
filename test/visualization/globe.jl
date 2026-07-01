using GeoMakie, Makie
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
            @test isa(sp, Makie.Surface)
            if colorbar
                @test isa(cb, Colorbar)
                @test cb.label[] == label
            else
                @test cb === nothing
            end
            @test ax.title[] == title
        end
    end

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
        @test isa(sp, Makie.Surface)
        @test cb === nothing
        @test_throws BoundsError globe_plot(FTLE, spectral_grid; time_index = 0)

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
        @test isa(sp, Makie.Surface)
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
        @test isa(sp, Makie.Surface)
        @test cb === nothing

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
        @test isa(sp, Makie.Surface)
        @test cb === nothing
    end
end
