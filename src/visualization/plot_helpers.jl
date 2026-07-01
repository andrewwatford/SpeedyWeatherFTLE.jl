function _plot_values(data)
    return data
end

_plot_values(result::FTLEResult) = result.ftle

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

export ftle_colorrange
