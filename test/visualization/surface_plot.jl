using GeoMakie, CairoMakie

@testset "surface_plot.jl" begin
    lons = -180:180
    lats = -90:90
    field = [exp(cosd(l)) + 3(y/90) for l in lons, y in lats]
    title = "Test Surface Plot"
    label = "Field Value"
    for coastlines in (true, false)
        for colorbar in (true, false)
            fig, ax, sp, cb = surface_plot(
                lons,
                lats,
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
