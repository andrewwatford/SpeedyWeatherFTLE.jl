using Documenter
using SpeedyWeatherFTLE

makedocs(
    sitename = "SpeedyWeatherFTLE",
    modules = [SpeedyWeatherFTLE],
    pages = [
        "Home" => "index.md",
    ],
)

deploydocs(
    repo = "github.com/andrewwatford/SpeedyWeatherFTLE.jl",
)