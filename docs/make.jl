using Documenter
using SpeedyWeatherFTLE

makedocs(
    sitename = "SpeedyWeatherFTLE",
    modules = [SpeedyWeatherFTLE],
    pages = [
        "Home" => "index.md",
        "Greeting" => "greet.md",
        "2D Turbulence" => "2d-turbulence.md",
        "Plotting the FTLE" => "get_FTLE.md",
    ],
)

deploydocs(
    repo = "github.com/andrewwatford/SpeedyWeatherFTLE.jl",
)