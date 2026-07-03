using Documenter
using SpeedyWeatherFTLE

makedocs(
    sitename = "SpeedyWeatherFTLE",
    modules = [SpeedyWeatherFTLE],
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
    ],
    checkdocs = :exports,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        edit_link = "main",
    ),
)
