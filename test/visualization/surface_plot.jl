using GeoMakie
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

    fig, ax, sp, cb = surface_plot(
        field;
        lon = -180:60:180,
        lat = -90:30:90,
        coastlines = false,
        colorbar = false,
    )

    @test isa(fig, Figure)
    @test isa(ax, GeoAxis)
    @test isa(sp, GeoMakie.Surface)
    @test cb === nothing

    @testset "FTLE matrix overload" begin
        spectral_grid = SpectralGrid(nlayers=1, trunc=6, Grid=FullGaussianGrid)
        FTLE = rand(spectral_grid.npoints, 4)
        FTLE_with_nan = copy(FTLE)
        FTLE_with_nan[1] = NaN
        field = ftle_field(FTLE, spectral_grid)
        final_field = ftle_field(FTLE[:, end], spectral_grid)
        result_kwargs = (;
            dist_km = 10,
            backwards = false,
            dynamics = false,
            rint_hours = 1,
        )
        result = FTLEResult(
            FTLE,
            spectral_grid,
            collect(0.0:3.0);
            result_kwargs...,
        )
        diagnostic_result = FTLEResult(
            FTLE,
            nothing,
            collect(0.0:3.0);
            result_kwargs...,
        )
        @test_throws DimensionMismatch FTLEResult(FTLE[:, end], spectral_grid, [3.0]; result_kwargs...)
        @test_throws DimensionMismatch FTLEResult(FTLE[1:end - 1, :], spectral_grid, collect(0.0:3.0); result_kwargs...)
        @test_throws DimensionMismatch FTLEResult(FTLE, spectral_grid, [0.0, 1.0]; result_kwargs...)
        @test_throws ArgumentError FTLEResult(FTLE, spectral_grid, [NaN, 1.0, 2.0, 3.0]; result_kwargs...)
        @test_throws ArgumentError FTLEResult(FTLE, spectral_grid, [0.0, Inf, 2.0, 3.0]; result_kwargs...)
        @test_throws ArgumentError FTLEResult(FTLE, spectral_grid, Any[0.0, "1.0", 2.0, 3.0]; result_kwargs...)
        @test_throws ArgumentError FTLEResult(
            FTLE,
            spectral_grid,
            collect(0.0:3.0);
            dist_km = 0,
            backwards = false,
            dynamics = false,
            rint_hours = 1,
        )
        @test_throws ArgumentError FTLEResult(
            FTLE,
            spectral_grid,
            collect(0.0:3.0);
            dist_km = 10,
            backwards = false,
            dynamics = false,
            rint_hours = Inf,
        )

        @test isa(field, Field)
        @test size(field) == size(FTLE)
        @test isa(final_field, Field)
        @test size(final_field) == (spectral_grid.npoints,)
        @test size(result) == size(FTLE)
        @test final_ftle(result) == FTLE[:, end]
        @test final_ftle_field(result) == final_field
        @test ftle_field(result; time_indices = 2) == ftle_field(view(FTLE, :, 2), spectral_grid)
        @test ftle_field(result; time_hour = 1.6) == ftle_field(view(FTLE, :, 3), spectral_grid)
        @test ftle_field(result; time_indices = :last) == final_field
        @test ftle_field(result; time_indices = :nonzero) == ftle_field(view(FTLE, :, 2:4), spectral_grid)
        @test_throws DimensionMismatch ftle_field(FTLE[1:end - 1, end], spectral_grid)
        @test_throws ArgumentError ftle_field(result; time_indices = :middle)
        @test_throws ArgumentError ftle_field(result; time_indices = 2, time_hour = 2.0)
        nonfinite_time_error = try
            ftle_field(result; time_hour = NaN)
            nothing
        catch err
            err
        end
        @test nonfinite_time_error isa ArgumentError
        @test occursin("time_hour must be finite", sprint(showerror, nonfinite_time_error))
        @test_throws ArgumentError ftle_field(diagnostic_result)
        @test_throws ArgumentError final_ftle_field(diagnostic_result)
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
        @test_throws ArgumentError surface_plot(FTLE, spectral_grid; time_index = 2, time_hour = 2.0, time_hours = result.time_hours)
        @test_throws ArgumentError surface_plot(FTLE, spectral_grid; time_hour = 2.0)
        @test_throws ArgumentError surface_plot(FTLE, spectral_grid; time_hour = NaN, time_hours = result.time_hours)
        @test_throws DimensionMismatch surface_plot(FTLE, spectral_grid; time_hour = 2.0, time_hours = [0.0, 1.0])

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
            FTLE,
            spectral_grid;
            time_hours = result.time_hours,
            time_hour = 2.6,
            colorbar = true,
            coastlines = false,
        )

        @test isa(fig, Figure)
        @test isa(ax, GeoAxis)
        @test isa(sp, GeoMakie.Surface)
        @test isa(cb, Colorbar)
        @test sp.colorrange[] ≈ collect(ftle_colorrange(view(FTLE, :, 4)))

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
        @test sp.colorrange[] ≈ collect(ftle_colorrange(FTLE[:, end]))

        fig, ax, sp, cb = surface_plot(
            FTLE_with_nan,
            spectral_grid;
            time_index = 1,
            colorbar = true,
            coastlines = false,
        )

        @test isa(fig, Figure)
        @test isa(ax, GeoAxis)
        @test isa(sp, GeoMakie.Surface)
        @test isa(cb, Colorbar)
        @test sp.colorrange[] ≈ collect(ftle_colorrange(view(FTLE_with_nan, :, 1)))

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
        @test sp.colorrange[] ≈ collect(ftle_colorrange(view(FTLE, :, 3)))
        @test_throws ArgumentError surface_plot(diagnostic_result)

        fig, ax, sp, cb = surface_plot(
            result;
            time_hour = 2.6,
            colorbar = true,
            coastlines = false,
        )

        @test isa(fig, Figure)
        @test isa(ax, GeoAxis)
        @test isa(sp, GeoMakie.Surface)
        @test isa(cb, Colorbar)
        @test sp.colorrange[] ≈ collect(ftle_colorrange(view(FTLE, :, 4)))
        @test_throws ArgumentError surface_plot(result; time_index = 2, time_hour = 2.0)
        @test_throws ArgumentError surface_plot(result; time_hour = Inf)

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
