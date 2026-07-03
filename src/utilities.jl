@inline _wrapped_lon_diff(lon1, lon2) = mod(lon1 - lon2 + 180, 360) - 180

function _check_dist_km(dist_km)
    dist_km isa Real && isfinite(dist_km) && dist_km > 0 ||
        throw(ArgumentError("dist_km must be finite and positive"))
    return Float64(dist_km)
end

function _check_radius(radius)
    radius isa Real && isfinite(radius) && radius > 0 ||
        throw(ArgumentError("radius must be finite and positive"))
    return Float64(radius)
end

function _check_positive(value, name)
    value isa Real && isfinite(value) && value > 0 ||
        throw(ArgumentError("$name must be finite and positive"))
    return Float64(value)
end

function _duration_magnitude(value, name)
    value isa Real && isfinite(value) || throw(ArgumentError("$name must be finite"))
    return abs(Float64(value))
end

function _check_time_hours(time_hours)
    for i in eachindex(time_hours)
        time_hours[i] isa Real && isfinite(time_hours[i]) ||
            throw(ArgumentError("time_hours[$i] must be real and finite"))
    end
    return nothing
end

_spatial_grid(grid) = grid
_spatial_grid(spectral_grid::SpectralGrid) = spectral_grid.grid

function _grid_npoints(npoints::Integer)
    npoints isa Bool && throw(ArgumentError("number of grid points must be an integer, not Bool"))
    npoints > 0 || throw(ArgumentError("number of grid points must be positive"))
    return Int(npoints)
end

_grid_npoints(grid) = length(first(RingGrids.get_londlatds(grid)))
_grid_npoints(spectral_grid::SpectralGrid) = spectral_grid.npoints

function _check_time_index(index::Integer, n_times)
    index isa Bool && throw(ArgumentError("time_indices must not contain Bool values"))
    1 <= index <= n_times || throw(BoundsError(1:n_times, index))
    return nothing
end

_checked_time_indices(time_indices, time_hours::AbstractVector) =
    _checked_time_indices(time_indices, length(time_hours), time_hours)

_checked_time_indices(::Colon, n_times, time_hours) = Base.OneTo(n_times)

function _checked_time_indices(selector::Symbol, n_times, time_hours)
    selector === :all && return Base.OneTo(n_times)
    selector === :first && return (1,)
    selector in (:last, :final) && return (n_times,)
    selector === :nonzero && return findall(t -> isfinite(t) && !iszero(t), time_hours)
    selector === :positive && return findall(t -> isfinite(t) && t > 0, time_hours)
    throw(ArgumentError("unsupported time_indices selector :$selector"))
end

function _checked_time_indices(index::Integer, n_times, time_hours)
    _check_time_index(index, n_times)
    return (index,)
end

_checked_time_indices(::Real, n_times, time_hours) =
    throw(ArgumentError("time_indices must be :, a Symbol, an integer, or integer indices"))

function _checked_time_indices(mask::AbstractVector{Bool}, n_times, time_hours)
    length(mask) == n_times ||
        throw(DimensionMismatch("boolean time_indices mask must have length $n_times"))
    return findall(mask)
end

function _checked_time_indices(indices, n_times, time_hours)
    for index in indices
        index isa Integer || throw(ArgumentError("time_indices must contain integer indices"))
        _check_time_index(index, n_times)
    end
    return indices
end

_selected_time_hours(time_hours, indices) = [Float64(time_hours[index]) for index in indices]

function _nearest_time_index(times, time_hour)
    isempty(times) && throw(ArgumentError("no time_hours are available"))
    isfinite(time_hour) || throw(ArgumentError("time_hour must be finite"))
    distances = abs.(times .- time_hour)
    return first(eachindex(times)) + argmin(distances) - 1
end

function _resolve_time_index(values; time_index, time_hour, time_hours)
    time_index !== nothing && time_hour !== nothing &&
        throw(ArgumentError("pass either time_index or time_hour, not both"))
    if time_hour !== nothing
        time_hours === nothing && throw(ArgumentError("time_hour requires time_hours"))
        length(time_hours) == size(values, 2) ||
            throw(DimensionMismatch("time_hours must contain one value per FTLE column"))
        time_index = _nearest_time_index(time_hours, time_hour)
    elseif time_index === nothing
        time_index = size(values, 2)
    end
    _check_time_index(time_index, size(values, 2))
    return time_index
end
