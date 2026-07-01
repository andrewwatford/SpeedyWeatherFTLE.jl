using GeoMakie, Makie
using SpeedyWeather
using SpeedyWeatherFTLE

@testset "surface_plot.jl" begin
    # Create a mock Field
    grid = HEALPixGrid(20)
    field = rand(grid)
    title = "Test Surface Plot"
    label = "Field Value"
    for coastlines in (true, false)
        for colorbar in (true, false)
            fig, ax, sp, cb = surface_plot(
                field;
                title = title,
                coastlines = coastlines,
                colorbar = colorbar,
                label = label,
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
        FTLE = rand(spectral_grid.npoints, 4)
        FTLE_with_nan = copy(FTLE)
        FTLE_with_nan[1] = NaN
        field = ftle_field(FTLE, spectral_grid)
        final_field = ftle_field(FTLE[:, end], spectral_grid)
        result = FTLEResult(
            FTLE,
            spectral_grid,
            collect(0.0:3.0);
            dist_km = 10,
            backwards = false,
            dynamics = false,
            rint_hours = 1,
        )

        @test isa(field, Field)
        @test size(field) == size(FTLE)
        @test isa(final_field, Field)
        @test size(final_field) == (spectral_grid.npoints,)
        @test size(result) == size(FTLE)
        @test final_ftle(result) == FTLE[:, end]
        @test final_ftle_field(result) == final_field
        @test ftle_field(result; time_indices = 2) == ftle_field(view(FTLE, :, 2), spectral_grid)
        @test ftle_field(result; time_indices = :last) == final_field
        @test ftle_field(result; time_indices = :nonzero) == ftle_field(view(FTLE, :, 2:4), spectral_grid)
        @test_throws DimensionMismatch ftle_field(FTLE[1:end - 1, end], spectral_grid)
        @test_throws ArgumentError ftle_field(result; time_indices = :middle)
        finite_FTLE = filter(isfinite, vec(FTLE_with_nan))
        @test ftle_colorrange(FTLE_with_nan) == (minimum(finite_FTLE), maximum(finite_FTLE))
        @test ftle_colorrange(result) == ftle_colorrange(FTLE)
        @test ftle_colorrange([-2.0, 1.0]; symmetric = true) == (-2.0, 2.0)
        @test collect(ftle_colorrange([1.0]; pad = 0.2)) ≈ [0.94, 1.06]
        @test_throws ArgumentError ftle_colorrange(FTLE; pad = -0.1)

        fig, ax, sp, cb = surface_plot(
            FTLE,
            spectral_grid;
            time_index = 4,
            title = title,
            colorbar = true,
            label = label,
            coastlines = false,
        )

        @test isa(fig, Figure)
        @test isa(ax, GeoAxis)
        @test isa(sp, GeoMakie.Surface)
        @test isa(cb, Colorbar)
        @test cb.label[] == label
        @test ax.title[] == title
        @test_throws BoundsError surface_plot(FTLE, spectral_grid; time_index = 0)

        shared_colorrange = (0.0, 1.0)
        fig, ax, sp, cb = surface_plot(
            FTLE,
            spectral_grid;
            time_index = 4,
            colorrange = shared_colorrange,
            coastlines = false,
            axis_kwargs = (; xlabel = "longitude"),
            surface_kwargs = (; transparency = false),
            colorbar_kwargs = (; vertical = true),
        )

        @test sp.colorrange[] == collect(shared_colorrange)
        @test ax.xlabel[] == "longitude"

        fig, ax, sp, cb = surface_plot(
            FTLE[:, end],
            spectral_grid;
            colorbar = true,
            coastlines = false,
        )

        @test isa(fig, Figure)
        @test isa(ax, GeoAxis)
        @test isa(sp, GeoMakie.Surface)
        @test isa(cb, Colorbar)
        @test cb.label[] == "FTLE [1/h]"

        fig, ax, sp, cb = surface_plot(
            FTLE[:, end],
            spectral_grid;
            title = title,
            colorbar = false,
            coastlines = false,
        )

        @test isa(fig, Figure)
        @test isa(ax, GeoAxis)
        @test isa(sp, GeoMakie.Surface)
        @test cb === nothing
        @test ax.title[] == title

        fig, ax, sp, cb = surface_plot(
            result;
            time_index = 3,
            colorbar = true,
            coastlines = false,
        )

        @test isa(fig, Figure)
        @test isa(ax, GeoAxis)
        @test isa(sp, GeoMakie.Surface)
        @test isa(cb, Colorbar)
        @test cb.label[] == "FTLE [1/h]"

        fig, ax, sp, cb = surface_plot(
            result;
            time_index = 3,
            colorbar = false,
            coastlines = false,
        )

        @test isa(fig, Figure)
        @test isa(ax, GeoAxis)
        @test isa(sp, GeoMakie.Surface)
        @test cb === nothing
    end
end
