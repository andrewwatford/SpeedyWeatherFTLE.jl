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
- `spectral_grid`: SpeedyWeather spectral grid used by the run, or `nothing`
  when only non-spatial diagnostics are needed.
- `time_hours`: selected tracker output times, measured in hours since release,
  with one entry per `ftle` column.
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
    ndims(ftle) == 2 ||
        throw(DimensionMismatch("ftle must be a matrix with dimensions (grid point, selected time)"))
    if spectral_grid !== nothing
        npoints = _grid_npoints(spectral_grid)
        size(ftle, 1) == npoints ||
            throw(DimensionMismatch("ftle has $(size(ftle, 1)) rows, but the spectral grid has $npoints grid points"))
    end
    length(time_hours) == size(ftle, 2) ||
        throw(DimensionMismatch("time_hours has length $(length(time_hours)), but ftle has $(size(ftle, 2)) time columns"))

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

function _nearest_time_index(times, time_hour::Real)
    isempty(times) && throw(ArgumentError("no time_hours are available"))
    isfinite(time_hour) ||
        throw(ArgumentError("time_hour must be finite, got $(time_hour)"))
    nearest_index = firstindex(times)
    nearest_distance = abs(times[nearest_index] - time_hour)
    for index in Iterators.drop(eachindex(times), 1)
        distance = abs(times[index] - time_hour)
        if distance < nearest_distance
            nearest_index = index
            nearest_distance = distance
        end
    end
    return nearest_index
end

function _require_spectral_grid(result::FTLEResult, operation)
    result.spectral_grid !== nothing && return result.spectral_grid
    throw(ArgumentError(
        "$(operation) requires result.spectral_grid, but this FTLEResult has no spectral grid. " *
        "Diagnostic-only FTLEResult objects cannot be converted to spatial fields or plotted; " *
        "pass FTLE data with an explicit grid, or use a result created with a spectral grid.",
    ))
end

function _resolve_time_index(FTLE_grid_time::AbstractMatrix; time_index, time_hour, time_hours)
    if time_index !== nothing && time_hour !== nothing
        throw(ArgumentError("pass either time_index or time_hour, not both"))
    end

    if time_hour !== nothing
        time_hours === nothing &&
            throw(ArgumentError("time_hour requires time_hours for matrix inputs"))
        length(time_hours) == size(FTLE_grid_time, 2) ||
            throw(DimensionMismatch("time_hours has length $(length(time_hours)), but FTLE_grid_time has $(size(FTLE_grid_time, 2)) time steps"))
        resolved_index = _nearest_time_index(time_hours, time_hour)
    elseif time_index === nothing
        resolved_index = size(FTLE_grid_time, 2)
    else
        resolved_index = time_index
    end

    1 <= resolved_index <= size(FTLE_grid_time, 2) ||
        throw(BoundsError(FTLE_grid_time, (:, resolved_index)))
    return resolved_index
end

function ftle_field(result::FTLEResult; time_indices=Colon(), time_hour::Union{Nothing,Real}=nothing)
    spectral_grid = _require_spectral_grid(result, "ftle_field")
    if time_hour !== nothing
        isequal(time_indices, Colon()) ||
            throw(ArgumentError("pass either time_indices or time_hour, not both"))
        selected_time_indices = [_resolve_time_index(
            result.ftle;
            time_index=nothing,
            time_hour,
            time_hours=result.time_hours,
        )]
    else
        selected_time_indices = _checked_time_indices(time_indices, result.time_hours)
    end

    if length(selected_time_indices) == 1
        return ftle_field(view(result.ftle, :, first(selected_time_indices)), spectral_grid)
    else
        return ftle_field(view(result.ftle, :, selected_time_indices), spectral_grid)
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
function final_ftle_field(result::FTLEResult)
    spectral_grid = _require_spectral_grid(result, "final_ftle_field")
    return ftle_field(final_ftle(result), spectral_grid)
end

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
