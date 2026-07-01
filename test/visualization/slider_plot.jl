using GeoMakie, Makie
using SpeedyWeather
using SpeedyWeatherFTLE

@testset "slider_plot.jl" begin
    # Create a mock Field
    grid = HEALPixGrid(20)
    field_ts = rand(grid, 5)  # A generic time-dependent field for plotting coverage.
    title = "Test Slider Plot"
    label = "Field Value"
    for coastlines in (true, false)
        for colorbar in (true, false)
            fig, ax, sp, cb = slider_plot(
                collect(1:5),
                field_ts,;
                title = title,
                coastlines = coastlines,
                coastline_color = :white,
                coastline_linewidth = 2,
                coastline_kwargs = (; overdraw = true),
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
            colorbar = true,
            coastlines = false,
        )

        @test isa(fig, Figure)
        @test isa(ax, GeoAxis)
        @test isa(sp, GeoMakie.Surface)
        @test isa(cb, Colorbar)
        @test cb.label[] == "FTLE [1/h]"

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

        handle = slider_plot(
            result;
            title = title,
            colorbar = false,
            coastlines = false,
            return_handle = true,
        )

        @test isa(handle, SliderPlotHandle)
        @test isa(handle.fig, Figure)
        @test isa(handle.ax, GeoAxis)
        @test isa(handle.sp, GeoMakie.Surface)
        @test handle.cb === nothing
        @test isa(handle.time_label, Label)
        @test handle.times == collect(1.0:4.0)
        @test handle.slider.value[] == 1
        @test set_slider_time!(handle, 2.6) === handle
        @test handle.slider.value[] == 3

        compact_handle = slider_plot(
            result;
            title = title,
            colorbar = false,
            coastlines = false,
            return_handle = true,
            time_label = false,
        )

        @test isa(compact_handle, SliderPlotHandle)
        @test compact_handle.time_label === nothing

        recorded_frames = Int[]
        fake_record(callback, fig, path, frames; framerate, kwargs...) = begin
            @test fig === handle.fig
            @test path == "ftle-slider.gif"
            @test framerate == 4
            for frame in frames
                push!(recorded_frames, frame)
                callback(frame)
            end
            path
        end

        returned_path = SpeedyWeatherFTLE._record_slider_animation!(
            fake_record,
            "ftle-slider.gif",
            handle,
            1:3;
            framerate = 4,
            record_kwargs = (;),
        )

        @test returned_path == "ftle-slider.gif"
        @test recorded_frames == [1, 2, 3]
        @test handle.slider.value[] == 3

        public_recorded_frames = Int[]
        public_fake_record(callback, fig, path, frames; framerate, kwargs...) = begin
            @test isa(fig, Figure)
            @test path == "public-ftle-slider.gif"
            @test framerate == 5
            for frame in frames
                push!(public_recorded_frames, frame)
                callback(frame)
            end
            path
        end

        public_returned_path = animate_slider_plot(
            "public-ftle-slider.gif",
            result;
            frames = 1:2,
            framerate = 5,
            record_function = public_fake_record,
            colorbar = false,
            coastlines = false,
        )

        @test public_returned_path == "public-ftle-slider.gif"
        @test public_recorded_frames == [1, 2]
    end
end
