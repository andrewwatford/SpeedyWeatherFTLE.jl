"""
    FTLEResult

Container returned by [`get_FTLE`](@ref), [`positive_FTLE`](@ref), or
[`negative_FTLE`](@ref) when `return_result = true`.

The result stores the FTLE matrix, the `SpectralGrid` used for the particle
tracking run, selected output times in hours, and metadata needed for plotting
or reusing saved particle files. `FTLEResult` supports array-like `size`,
`axes`, `length`, `eachindex`, `eltype`, and `getindex` by forwarding to its
`ftle` matrix.

# Fields

- `ftle`: matrix with dimensions `(grid point, selected time)`.
- `spectral_grid`: SpeedyWeather spectral grid used by the run.
- `time_hours`: selected tracker output times, measured in hours since release.
- `particle_file_path`: path to the saved particle file, or `nothing`.
- `dist_km`: initial particle perturbation distance in kilometres.
- `backwards`: whether the simulation ran backward in time.
- `direction`: `:positive` or `:negative`.
- `dynamics`: whether SpeedyWeather dynamics were enabled.
- `rint_hours`: particle-tracker output cadence in hours.
"""
struct FTLEResult{F, S, T, P}
    ftle::F
    spectral_grid::S
    time_hours::T
    particle_file_path::P
    dist_km::Float64
    backwards::Bool
    direction::Symbol
    dynamics::Bool
    rint_hours::Float64
end

function FTLEResult(
    ftle,
    spectral_grid,
    time_hours;
    particle_file_path=nothing,
    dist_km,
    backwards,
    direction=backwards ? :negative : :positive,
    dynamics,
    rint_hours,
)
    return FTLEResult(
        ftle,
        spectral_grid,
        time_hours,
        particle_file_path,
        Float64(dist_km),
        Bool(backwards),
        Symbol(direction),
        Bool(dynamics),
        Float64(rint_hours),
    )
end

function ftle_field(result::FTLEResult; time_indices=Colon())
    selected_time_indices = _checked_time_indices(time_indices, result.time_hours)
    if length(selected_time_indices) == 1
        return ftle_field(view(result.ftle, :, first(selected_time_indices)), result.spectral_grid)
    else
        return ftle_field(view(result.ftle, :, selected_time_indices), result.spectral_grid)
    end
end

"""
    final_ftle(result::FTLEResult)

Return the final selected FTLE column from `result.ftle` as a vector.
"""
final_ftle(result::FTLEResult) = result.ftle[:, end]

"""
    final_ftle_field(result::FTLEResult)

Return the final selected FTLE column as a `RingGrids.Field`, ready for
interpolation or plotting.
"""
final_ftle_field(result::FTLEResult) = ftle_field(result; time_indices=:last)

Base.size(result::FTLEResult) = size(result.ftle)
Base.size(result::FTLEResult, dim::Integer) = size(result.ftle, dim)
Base.axes(result::FTLEResult) = axes(result.ftle)
Base.axes(result::FTLEResult, dim::Integer) = axes(result.ftle, dim)
Base.length(result::FTLEResult) = length(result.ftle)
Base.eachindex(result::FTLEResult) = eachindex(result.ftle)
Base.eltype(result::Type{<:FTLEResult{F}}) where {F} = eltype(F)
Base.getindex(result::FTLEResult, indices...) = getindex(result.ftle, indices...)

function _time_hours_summary(time_hours)
    if isempty(time_hours)
        return "none"
    elseif length(time_hours) == 1
        return "$(only(time_hours)) h"
    else
        return "$(first(time_hours)) to $(last(time_hours)) h"
    end
end

function Base.show(io::IO, result::FTLEResult)
    npoints, ntimes = size(result.ftle)
    print(
        io,
        "FTLEResult(",
        npoints,
        " grid points, ",
        ntimes,
        " times, direction=",
        result.direction,
        ", time_hours=",
        _time_hours_summary(result.time_hours),
        ")",
    )
end

function Base.show(io::IO, ::MIME"text/plain", result::FTLEResult)
    npoints, ntimes = size(result.ftle)
    println(io, "FTLEResult")
    println(io, "  grid points: ", npoints)
    println(io, "  times: ", ntimes, " (", _time_hours_summary(result.time_hours), ")")
    println(io, "  direction: ", result.direction)
    println(io, "  dist_km: ", result.dist_km)
    println(io, "  dynamics: ", result.dynamics)
    println(io, "  rint_hours: ", result.rint_hours)
    print(io, "  particle_file_path: ", result.particle_file_path === nothing ? "nothing" : result.particle_file_path)
end

export FTLEResult
export final_ftle
export final_ftle_field
