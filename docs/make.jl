using Documenter
using SpeedyWeatherFTLE

makedocs(
    sitename = "SpeedyWeatherFTLE",
    modules = [SpeedyWeatherFTLE],
    pages = [
        "Home" => "index.md",
        "Greeting" => "greet.md",
        "2D Turbulence" => "2d-turbulence.md",
        "Lyapunov Exponents and FTLE" => "lyapunov-FTLE.md",
    ],
)

deploydocs(
    repo = "github.com/andrewwatford/SpeedyWeatherFTLE.jl",
)