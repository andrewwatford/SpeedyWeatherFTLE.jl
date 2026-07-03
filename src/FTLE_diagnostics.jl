"""
    stretching_factor!(stretch, ftle, time_hours)

In-place conversion from FTLE values to finite-time stretching factors. The
time argument must use the same unit as the FTLE rate.
"""
function stretching_factor!(stretch, ftle::AbstractVector, time_hour::Real)
    length(stretch) == length(ftle) ||
        throw(DimensionMismatch("stretch must have length $(length(ftle))"))
    duration = _duration_magnitude(time_hour, "time_hour")
    @inbounds for i in eachindex(stretch, ftle)
        stretch[i] = exp(ftle[i] * duration)
    end
    return stretch
end

"""
    stretching_factor(ftle, time_hours)
    stretching_factor(result::FTLEResult)

Convert FTLE values to finite-time stretching factors with
`exp.(ftle .* abs.(time))`.
"""
function stretching_factor(ftle::AbstractVector, time_hour::Real)
    return stretching_factor!(similar(ftle, Float64), ftle, time_hour)
end

function stretching_factor!(stretch, ftle::AbstractMatrix, time_hours::AbstractVector{<:Real})
    size(stretch) == size(ftle) ||
        throw(DimensionMismatch("stretch must have size $(size(ftle))"))
    length(time_hours) == size(ftle, 2) ||
        throw(DimensionMismatch("time_hours must contain one value per FTLE column"))

    @inbounds for (j, time_hour) in enumerate(time_hours)
        duration = _duration_magnitude(time_hour, "time_hours[$j]")
        for i in axes(ftle, 1)
            stretch[i, j] = exp(ftle[i, j] * duration)
        end
    end
    return stretch
end

function stretching_factor(ftle::AbstractMatrix, time_hours::AbstractVector{<:Real})
    return stretching_factor!(similar(ftle, Float64), ftle, time_hours)
end

stretching_factor(result::FTLEResult) = stretching_factor(result.ftle, result.time_hours)
