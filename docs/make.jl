using CairoMakie
using Documenter
using SpeedyWeatherFTLE

makedocs(
    sitename = "SpeedyWeatherFTLE",
    modules = [SpeedyWeatherFTLE],
    pages = [
        "Home" => "index.md",
        "Concepts and Data Layout" => "concepts.md",
        "Running Simulations" => "simulation.md",
        "Particle Files" => "particle_files.md",
        "Plotting" => "plotting.md",
        "API Reference" => "api.md",
    ],
    checkdocs = :exports,
    format = Documenter.HTML(prettyurls = true),
)

deploydocs(
    repo = "github.com/andrewwatford/SpeedyWeatherFTLE.jl",
)
