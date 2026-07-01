"""
    stretching_factor(ftle, time_hours)
    stretching_factor(result::FTLEResult)

Convert FTLE values to finite-time stretching factors.

For a scalar integration time `T`, this returns `exp.(ftle .* T)`. For an
`FTLE_grid_time` matrix, `time_hours` must contain one time per column and the
conversion is applied column by column. Values are dimensionless; `time_hours`
must use the same time unit used to compute `ftle`, which is hours for
SpeedyWeatherFTLE particle post-processing.

`stretching_factor(result)` uses `result.ftle` and `result.time_hours`.
"""
function stretching_factor end

"""
    stretching_factor!(stretch, ftle, time_hours)

In-place form of [`stretching_factor`](@ref), writing into `stretch` and
returning it.
"""
function stretching_factor! end

function stretching_factor!(stretch, ftle::AbstractVector, time_hour::Real)
    length(stretch) == length(ftle) ||
        throw(DimensionMismatch("stretch must have length $(length(ftle))"))

    @inbounds for i in eachindex(stretch, ftle)
        stretch[i] = exp(ftle[i] * time_hour)
    end

    return stretch
end

function stretching_factor(ftle::AbstractVector, time_hour::Real)
    stretch = similar(ftle, Float64)
    return stretching_factor!(stretch, ftle, time_hour)
end

function stretching_factor!(stretch, ftle::AbstractMatrix, time_hours::AbstractVector{<:Real})
    size(stretch) == size(ftle) ||
        throw(DimensionMismatch("stretch must have size $(size(ftle))"))
    size(ftle, 2) == length(time_hours) ||
        throw(DimensionMismatch("time_hours must contain one time per FTLE column"))

    col_axis = axes(ftle, 2)
    @inbounds for (time_index, time_hour) in enumerate(time_hours)
        j = col_axis[time_index]
        for i in axes(ftle, 1)
            stretch[i, j] = exp(ftle[i, j] * time_hour)
        end
    end

    return stretch
end

function stretching_factor(ftle::AbstractMatrix, time_hours::AbstractVector{<:Real})
    stretch = similar(ftle, Float64)
    return stretching_factor!(stretch, ftle, time_hours)
end

stretching_factor(result::FTLEResult) = stretching_factor(result.ftle, result.time_hours)

export stretching_factor
export stretching_factor!
