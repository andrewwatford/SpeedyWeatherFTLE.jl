using RingGrids: Field, interpolate

const GeoMakie = Base.loaded_modules[Base.PkgId(Base.UUID("db073c08-6b98-4ee5-b6a4-5efafb3259c6"), "GeoMakie")]
const Makie = GeoMakie.Makie
const NoShading = Makie.NoShading
const Figure = Makie.Figure
const GeoAxis = GeoMakie.GeoAxis
const Colorbar = Makie.Colorbar
const Relative = Makie.Relative
const GridLayout = Makie.GridLayout
const SliderGrid = Makie.SliderGrid
const Label = Makie.Label
const lift = Makie.lift
const surface! = Makie.surface!
const lines! = Makie.lines!

include("surface_plot.jl")
include("slider_plot.jl")
include("globe.jl")
