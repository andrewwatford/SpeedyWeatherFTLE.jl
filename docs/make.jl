using Documenter
using SpeedyWeatherFTLE

makedocs(
    sitename = "SpeedyWeatherFTLE",
    modules = [SpeedyWeatherFTLE],
    pages = [
        "Home" => "index.md",
        "Greeting" => "greet.md",
    ],
)

deploydocs(
    repo = "github.com/andrewwatford/SpeedyWeatherFTLE.jl",
)