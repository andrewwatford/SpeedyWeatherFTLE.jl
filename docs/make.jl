using CairoMakie
using Documenter
using GeoMakie
using SpeedyWeatherFTLE

CairoMakie.activate!()

makedocs(
    sitename = "SpeedyWeatherFTLE",
    modules = [SpeedyWeatherFTLE],
    pages = [
        "Home" => "index.md",
        "Running FTLE" => "simulation.md",
        "Worked Example" => "zonal_jet.md",
        "Plotting" => "plotting.md",
        "API" => "api.md",
    ],
    checkdocs = :exports,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        edit_link = "main",
    ),
)
