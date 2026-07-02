using RingGrids
using SpeedyWeather
using SpeedyWeatherFTLE
using Test

@testset "FTLE_result.jl" begin
    spectral_grid = SpectralGrid(nlayers=1, trunc=6, Grid=FullGaussianGrid)
    FTLE = rand(spectral_grid.npoints, 4)
    FTLE_with_nan = copy(FTLE)
    FTLE_with_nan[1] = NaN
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

    field = ftle_field(FTLE, spectral_grid)
    final_field = ftle_field(FTLE[:, end], spectral_grid)

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
end
