"""
    FTLEResult(ftle, spectral_grid, time_hours; dist_km, backwards, ...)

Small named container returned by `FTLE(...; return_result = true)`.
"""
struct FTLEResult{F, S, T}
    ftle::F
    spectral_grid::S
    time_hours::T
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
    dist_km,
    backwards=false,
    direction=backwards ? :backward : :forward,
    dynamics=false,
    rint_hours=NaN,
)
    ndims(ftle) == 2 ||
        throw(DimensionMismatch("ftle must be a matrix with dimensions (grid point, time)"))
    spectral_grid === nothing || size(ftle, 1) == _grid_npoints(spectral_grid) ||
        throw(DimensionMismatch("ftle row count does not match spectral_grid"))
    length(time_hours) == size(ftle, 2) ||
        throw(DimensionMismatch("time_hours must contain one value per FTLE column"))
    _check_time_hours(time_hours)
    dist_km = _check_dist_km(dist_km)
    !isnan(rint_hours) && _check_positive(rint_hours, "rint_hours")

    return FTLEResult(
        ftle,
        spectral_grid,
        time_hours,
        dist_km,
        Bool(backwards),
        Symbol(direction),
        Bool(dynamics),
        Float64(rint_hours),
    )
end

Base.size(result::FTLEResult) = size(result.ftle)
Base.size(result::FTLEResult, dim::Integer) = size(result.ftle, dim)
Base.axes(result::FTLEResult) = axes(result.ftle)
Base.axes(result::FTLEResult, dim::Integer) = axes(result.ftle, dim)
Base.length(result::FTLEResult) = length(result.ftle)
Base.eachindex(result::FTLEResult) = eachindex(result.ftle)
Base.eltype(result::FTLEResult) = eltype(result.ftle)
Base.eltype(::Type{<:FTLEResult{F}}) where {F} = eltype(F)
Base.getindex(result::FTLEResult, indices...) = getindex(result.ftle, indices...)

function _require_spectral_grid(result::FTLEResult, operation)
    result.spectral_grid !== nothing && return result.spectral_grid
    throw(ArgumentError("$operation requires result.spectral_grid"))
end

"""
    final_ftle(result::FTLEResult)

Return the final selected FTLE column from `result`.
"""
final_ftle(result::FTLEResult) = result.ftle[:, end]

"""
    ftle_field(ftle, grid_or_spectral_grid)
    ftle_field(result::FTLEResult; time_indices = :, time_hour = nothing)

Convert FTLE vectors, matrices, or selected `FTLEResult` columns to
`RingGrids.Field` values.
"""
function ftle_field(ftle::AbstractMatrix, grid_or_spectral_grid)
    npoints = _grid_npoints(grid_or_spectral_grid)
    size(ftle, 1) == npoints ||
        throw(DimensionMismatch("ftle has $(size(ftle, 1)) rows, expected $npoints"))
    return Field(ftle, _spatial_grid(grid_or_spectral_grid))
end

function ftle_field(ftle::AbstractVector, grid_or_spectral_grid)
    npoints = _grid_npoints(grid_or_spectral_grid)
    length(ftle) == npoints ||
        throw(DimensionMismatch("ftle has length $(length(ftle)), expected $npoints"))
    return Field(ftle, _spatial_grid(grid_or_spectral_grid))
end

function ftle_field(result::FTLEResult; time_indices=Colon(), time_hour=nothing)
    spectral_grid = _require_spectral_grid(result, "ftle_field")
    if time_hour !== nothing
        isequal(time_indices, Colon()) ||
            throw(ArgumentError("pass either time_indices or time_hour, not both"))
        index = _resolve_time_index(result.ftle; time_index=nothing, time_hour, time_hours=result.time_hours)
        return ftle_field(view(result.ftle, :, index), spectral_grid)
    end

    selected = _checked_time_indices(time_indices, result.time_hours)
    length(selected) == 1 && return ftle_field(view(result.ftle, :, first(selected)), spectral_grid)
    return ftle_field(view(result.ftle, :, selected), spectral_grid)
end

"""
    final_ftle_field(result::FTLEResult)

Return the final selected FTLE column as a `RingGrids.Field`.
"""
final_ftle_field(result::FTLEResult) = ftle_field(final_ftle(result), _require_spectral_grid(result, "final_ftle_field"))

function Base.show(io::IO, result::FTLEResult)
    print(io, "FTLEResult(", size(result.ftle, 1), " points, ", size(result.ftle, 2), " times, direction=", result.direction, ")")
end
