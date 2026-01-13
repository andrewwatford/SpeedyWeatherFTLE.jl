using GeoMakie, CairoMakie
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
end