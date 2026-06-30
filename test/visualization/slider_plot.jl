using GeoMakie, Makie
using SpeedyWeather
using SpeedyWeatherFTLE

@testset "slider_plot.jl" begin
    # Create a mock Field
    grid = HEALPixGrid(20)
    field_ts = rand(grid, 5)  # Assuming a time series with 5 time steps
    title = "Test Slider Plot"
    label = "Field Value"
    for coastlines in (true, false)
        for colorbar in (true, false)
            fig, ax, sp, cb = slider_plot(
                collect(1:5),
                field_ts,;
                title = title,
                coastlines = coastlines,
                colorbar = colorbar,
                colorbar_label = label,
            )
            @test isa(fig, Figure)
            @test isa(ax, GeoAxis)
            @test isa(sp, GeoMakie.Surface)
            if colorbar
                @test isa(cb, Colorbar)
                @test cb.label[] == label
            end
            @test ax.title[] == title
        end
    end

    @testset "FTLE matrix overload" begin
        spectral_grid = SpectralGrid(nlayers=1, trunc=6, Grid=FullGaussianGrid)
        FTLE = rand(spectral_grid.npoints, 5)
        FTLE[:, 1] .= NaN
        result = FTLEResult(
            FTLE,
            spectral_grid,
            collect(0.0:4.0);
            dist_km = 10,
            backwards = false,
            dynamics = false,
            rint_hours = 1,
        )

        fig, ax, sp, cb = slider_plot(
            collect(0:4),
            FTLE,
            spectral_grid;
            title = title,
            colorbar = true,
            colorbar_label = label,
            coastlines = false,
        )

        @test isa(fig, Figure)
        @test isa(ax, GeoAxis)
        @test isa(sp, GeoMakie.Surface)
        @test isa(cb, Colorbar)
        @test cb.label[] == label
        @test ax.title[] == title
        @test_throws DimensionMismatch slider_plot(collect(0:3), FTLE, spectral_grid)
        @test_throws BoundsError slider_plot(collect(0:4), FTLE, spectral_grid; start_index = 0)

        fig, ax, sp, cb = slider_plot(
            collect(0:4),
            FTLE,
            spectral_grid;
            start_index = 3,
            title = title,
            colorbar = false,
            coastlines = false,
        )

        @test isa(fig, Figure)
        @test isa(ax, GeoAxis)
        @test isa(sp, GeoMakie.Surface)
        @test cb === nothing
        @test ax.title[] == title

        fig, ax, sp, cb = slider_plot(
            result;
            title = title,
            colorbar = false,
            coastlines = false,
        )

        @test isa(fig, Figure)
        @test isa(ax, GeoAxis)
        @test isa(sp, GeoMakie.Surface)
        @test cb === nothing
    end
end
