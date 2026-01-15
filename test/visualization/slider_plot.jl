using GeoMakie, GLMakie
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
                Vector{Float64}(1:5),
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
end