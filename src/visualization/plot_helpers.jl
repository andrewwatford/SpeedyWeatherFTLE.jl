function _plot_values(data)
    return data
end

_plot_values(result::FTLEResult) = result.ftle

const _FTLE_COLORBAR_LABEL = "FTLE [1/h]"
const _FTLE_SLIDER_LABEL = "Integration time [h]"

"""
    ftle_colorrange(data...; symmetric = false, pad = 0)

Return finite color limits spanning one or more FTLE arrays, fields, or
[`FTLEResult`](@ref) objects.

Use this when comparing multiple plots so their colors share the same scale.
With `symmetric = true`, the returned range is centred on zero. `pad` expands
the range by a fraction of its width, which can make nearly uniform fields more
legible.
"""
function ftle_colorrange(data...; symmetric::Bool=false, pad::Real=0)
    pad >= 0 || throw(ArgumentError("pad must be non-negative"))

    min_val = Inf
    max_val = -Inf
    found_finite = false
    for dataset in data
        for value in _plot_values(dataset)
            if isfinite(value)
                finite_value = float(value)
                min_val = min(min_val, finite_value)
                max_val = max(max_val, finite_value)
                found_finite = true
            end
        end
    end

    if !found_finite
        min_val, max_val = 0.0, 1.0
    elseif symmetric
        max_abs = max(abs(min_val), abs(max_val), eps(Float64))
        min_val, max_val = -max_abs, max_abs
    elseif min_val == max_val
        half_width = 0.05 * max(abs(min_val), 1.0)
        min_val, max_val = min_val - half_width, max_val + half_width
    end

    if pad > 0
        half_pad = pad * (max_val - min_val) / 2
        min_val -= half_pad
        max_val += half_pad
    end

    return (min_val, max_val)
end

function _resolve_colorrange(data, colorrange)
    if colorrange === nothing || colorrange === :auto
        return ftle_colorrange(data)
    elseif colorrange === :symmetric
        return ftle_colorrange(data; symmetric=true)
    else
        return colorrange
    end
end

function _nearest_time_index(times, time_hour::Real)
    isempty(times) && throw(ArgumentError("no time_hours are available"))
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

export ftle_colorrange
