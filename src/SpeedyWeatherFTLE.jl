module SpeedyWeatherFTLE

using Dates: Hour, Second
using NCDatasets
using RingGrids
using SpeedyWeather

export FTLEResult
export FTLE_from_particle_file, FTLE_from_particle_file!
export FTLE_from_particles, FTLE_from_particles!
export Re
export final_ftle, final_ftle_field
export ftle_field
export get_FTLE, positive_FTLE, negative_FTLE
export initial_FTLE_particle_positions, initial_FTLE_particle_positions!
export stretching_factor, stretching_factor!

"""
    Re

Average Earth radius used by SpeedyWeatherFTLE, in metres.
"""
const Re = 6.371e6
const _MAX_STENCIL_DEGREES = 20.0
const _PARTICLE_FILE_METADATA_PREFIX = "SpeedyWeatherFTLE_"
const _DIST_KM_ATTRIBUTE = _PARTICLE_FILE_METADATA_PREFIX * "dist_km"
const _PARTICLE_ORDER_ATTRIBUTE = _PARTICLE_FILE_METADATA_PREFIX * "particle_order"
const _PARTICLE_ORDER = "east,west,north,south"

@inline _wrapped_lon_diff(lon1, lon2) = mod(lon1 - lon2 + 180, 360) - 180

function _check_dist_km(dist_km)
    dist_km isa Real && isfinite(dist_km) && dist_km > 0 ||
        throw(ArgumentError("dist_km must be finite and positive"))
    return Float64(dist_km)
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

function _check_stencil_centers(londs, latds, dist_km)
    length(londs) == length(latds) ||
        throw(DimensionMismatch("londs and latds must have the same length"))

    dist_km = _check_dist_km(dist_km)
    delta_lat = rad2deg(dist_km * 1000 / Re)
    delta_lat <= _MAX_STENCIL_DEGREES ||
        throw(ArgumentError("dist_km=$dist_km is too large for the local spherical stencil"))

    for i in eachindex(londs, latds)
        lon = londs[i]
        lat = latds[i]
        lon isa Real && lat isa Real && isfinite(lon) && isfinite(lat) ||
            throw(ArgumentError("stencil centers must be finite real values"))
        abs(lat) + delta_lat < 90 ||
            throw(ArgumentError("dist_km=$dist_km makes the stencil cross a pole at latitude $lat"))
        delta_lon = delta_lat / abs(cosd(lat))
        isfinite(delta_lon) && delta_lon <= _MAX_STENCIL_DEGREES ||
            throw(ArgumentError("dist_km=$dist_km creates a nonlocal east/west stencil at latitude $lat"))
    end

    return length(londs), delta_lat
end

"""
    initial_FTLE_particle_positions!(plonds, platds, londs, latds, dist_km)
    initial_FTLE_particle_positions!(plonds, platds, grid_or_spectral_grid, dist_km)

Fill longitude and latitude arrays with the four-particle FTLE release stencil
around each grid point. Particles are ordered east, west, north, south.
"""
function initial_FTLE_particle_positions!(plonds, platds, londs, latds, dist_km)
    npoints, delta_lat = _check_stencil_centers(londs, latds, dist_km)
    expected_length = 4 * npoints
    length(plonds) == expected_length || throw(DimensionMismatch("plonds must have length $expected_length"))
    length(platds) == expected_length || throw(DimensionMismatch("platds must have length $expected_length"))

    @inbounds for i in 1:npoints
        p = 4 * i
        delta_lon = delta_lat / cosd(latds[i])
        plonds[p - 3] = londs[i] + delta_lon
        platds[p - 3] = latds[i]
        plonds[p - 2] = londs[i] - delta_lon
        platds[p - 2] = latds[i]
        plonds[p - 1] = londs[i]
        platds[p - 1] = latds[i] + delta_lat
        plonds[p] = londs[i]
        platds[p] = latds[i] - delta_lat
    end

    return plonds, platds
end

function initial_FTLE_particle_positions!(plonds, platds, grid_or_spectral_grid, dist_km)
    londs, latds = RingGrids.get_londlatds(_spatial_grid(grid_or_spectral_grid))
    return initial_FTLE_particle_positions!(plonds, platds, londs, latds, dist_km)
end

"""
    initial_FTLE_particle_positions(londs, latds, dist_km)
    initial_FTLE_particle_positions(grid_or_spectral_grid, dist_km)

Allocate and return longitude and latitude vectors for the canonical
east/west/north/south FTLE release stencil.
"""
function initial_FTLE_particle_positions(londs, latds, dist_km)
    plonds = Vector{Float64}(undef, 4 * length(londs))
    platds = similar(plonds)
    return initial_FTLE_particle_positions!(plonds, platds, londs, latds, dist_km)
end

function initial_FTLE_particle_positions(grid_or_spectral_grid, dist_km)
    npoints = _grid_npoints(grid_or_spectral_grid)
    plonds = Vector{Float64}(undef, 4 * npoints)
    platds = similar(plonds)
    return initial_FTLE_particle_positions!(plonds, platds, grid_or_spectral_grid, dist_km)
end

function perturb_positions_FTLE(particles, londs, latds, dist_km)
    npoints, delta_lat = _check_stencil_centers(londs, latds, dist_km)
    expected_length = 4 * npoints
    length(particles) == expected_length ||
        throw(DimensionMismatch("particles must have length $expected_length"))

    @inbounds for i in 1:npoints
        p = 4 * i
        delta_lon = delta_lat / cosd(latds[i])
        particles[p - 3] = Particle(londs[i] + delta_lon, latds[i])
        particles[p - 2] = Particle(londs[i] - delta_lon, latds[i])
        particles[p - 1] = Particle(londs[i], latds[i] + delta_lat)
        particles[p] = Particle(londs[i], latds[i] - delta_lat)
    end

    return particles
end

function displacement_gradient_matrix_central!(B, plonds, platds, dist_km)
    _check_dist_km(dist_km)
    length(plonds) == length(platds) ||
        throw(DimensionMismatch("plonds and platds must have the same length"))
    length(plonds) % 4 == 0 ||
        throw(ArgumentError("particle vectors must contain four particles per grid point"))

    npoints = length(plonds) ÷ 4
    size(B) == (2, 2, npoints) ||
        throw(DimensionMismatch("B must have size (2, 2, $npoints)"))

    dfac = Re / (1000 * dist_km) / 2
    @inbounds for i in 1:npoints
        p = 4 * i
        B[1, 1, i] = deg2rad(_wrapped_lon_diff(plonds[p - 3], plonds[p - 2])) *
                     cosd((platds[p - 3] + platds[p - 2]) / 2) * dfac
        B[2, 1, i] = deg2rad(platds[p - 3] - platds[p - 2]) * dfac
        B[1, 2, i] = deg2rad(_wrapped_lon_diff(plonds[p - 1], plonds[p])) *
                     cosd((platds[p - 1] + platds[p]) / 2) * dfac
        B[2, 2, i] = deg2rad(platds[p - 1] - platds[p]) * dfac
    end

    return B
end

function displacement_gradient_matrix_central(plonds, platds, dist_km)
    B = Array{Float64}(undef, 2, 2, length(plonds) ÷ 4)
    return displacement_gradient_matrix_central!(B, plonds, platds, dist_km)
end

@inline function _largest_cauchy_green_eigenvalue(a, b, c, d)
    c11 = a*a + c*c
    c12 = a*b + c*d
    c22 = b*b + d*d
    return (c11 + c22 + sqrt((c11 - c22)^2 + 4 * c12^2)) / 2
end

function FTLE_over_grid!(ftle, B, time)
    npoints = size(B, 3)
    length(ftle) == npoints || throw(DimensionMismatch("ftle must have length $npoints"))

    duration = _duration_magnitude(time, "time")
    if iszero(duration)
        fill!(ftle, NaN)
        return ftle
    end

    @inbounds for i in 1:npoints
        lmax = _largest_cauchy_green_eigenvalue(B[1, 1, i], B[1, 2, i], B[2, 1, i], B[2, 2, i])
        ftle[i] = log(lmax) / (2 * duration)
    end

    return ftle
end

function FTLE_over_grid(B, time)
    ftle = Vector{Float64}(undef, size(B, 3))
    return FTLE_over_grid!(ftle, B, time)
end

_particle_column(A::Union{NCDatasets.Variable, NCDatasets.CFVariable}, index) = A[:, index]
_particle_column(A::AbstractArray, index) = view(A, :, index)
_particle_column(A, index) = A[:, index]

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

function _check_particle_column(plonds, platds, time_index)
    length(plonds) == length(platds) ||
        throw(DimensionMismatch("lon and lat columns at time index $time_index have different lengths"))

    for i in eachindex(plonds, platds)
        lon = plonds[i]
        lat = platds[i]
        (ismissing(lon) || ismissing(lat)) &&
            throw(ArgumentError("particle positions at time index $time_index contain missing values"))
        lon isa Real && lat isa Real && isfinite(lon) && isfinite(lat) ||
            throw(ArgumentError("particle positions at time index $time_index must be finite real values"))
    end

    return nothing
end

"""
    FTLE_from_particles!(ftle, B, plonds_time, platds_time, time_hours, grid_or_npoints, dist_km; time_indices = :)

Compute FTLE from particle trajectories into caller-provided output and work
arrays. Longitude and latitude arrays must be `(particle, time)` with four
particles per grid point in east, west, north, south order.
"""
function FTLE_from_particles!(
    ftle,
    B,
    plonds_time,
    platds_time,
    time_hours::AbstractVector{<:Real},
    grid_or_npoints,
    dist_km;
    time_indices=Colon(),
)
    _check_dist_km(dist_km)
    size(plonds_time) == size(platds_time) ||
        throw(DimensionMismatch("plonds_time and platds_time must have the same size"))

    npoints = _grid_npoints(grid_or_npoints)
    expected_particles = 4 * npoints
    size(plonds_time, 1) == expected_particles ||
        throw(DimensionMismatch("particle arrays must have $expected_particles rows"))
    size(plonds_time, 2) == length(time_hours) ||
        throw(DimensionMismatch("particle arrays and time_hours must have the same number of time samples"))

    _check_time_hours(time_hours)
    selected = _checked_time_indices(time_indices, time_hours)
    size(ftle) == (npoints, length(selected)) ||
        throw(DimensionMismatch("ftle must have size ($npoints, $(length(selected)))"))
    size(B) == (2, 2, npoints) ||
        throw(DimensionMismatch("B must have size (2, 2, $npoints)"))

    for (out_index, time_index) in enumerate(selected)
        plonds = _particle_column(plonds_time, time_index)
        platds = _particle_column(platds_time, time_index)
        _check_particle_column(plonds, platds, time_index)
        displacement_gradient_matrix_central!(B, plonds, platds, dist_km)
        FTLE_over_grid!(view(ftle, :, out_index), B, time_hours[time_index])
    end

    return ftle
end

"""
    FTLE_from_particles(plonds_time, platds_time, time_hours, grid_or_npoints, dist_km; time_indices = :)

Allocate and compute FTLE from in-memory particle trajectories. Returns
`ftle, selected_time_hours`.
"""
function FTLE_from_particles(plonds_time, platds_time, time_hours::AbstractVector{<:Real}, grid_or_npoints, dist_km; time_indices=Colon())
    npoints = _grid_npoints(grid_or_npoints)
    _check_time_hours(time_hours)
    selected = _checked_time_indices(time_indices, time_hours)
    ftle = Array{Float64}(undef, npoints, length(selected))
    B = Array{Float64}(undef, 2, 2, npoints)
    FTLE_from_particles!(ftle, B, plonds_time, platds_time, time_hours, npoints, dist_km; time_indices=selected)
    return ftle, _selected_time_hours(time_hours, selected)
end

function _particle_file_time_hours(ds)
    times = ds["time"][:]
    isempty(times) && throw(ArgumentError("particle file time variable is empty"))

    first_time = first(times)
    if first_time isa Real
        _check_time_hours(times)
        return Float64.(times .- first_time)
    end

    try
        return Float64.((times .- first_time) ./ Hour(1))
    catch
        throw(ArgumentError("particle file time variable must contain real elapsed times or date/time values"))
    end
end

function _validate_particle_file_metadata(ds, dist_km)
    if haskey(ds.attrib, _DIST_KM_ATTRIBUTE)
        file_dist_km = ds.attrib[_DIST_KM_ATTRIBUTE]
        file_dist_km isa Real ||
            throw(ArgumentError("particle file dist_km metadata must be real-valued"))
        isapprox(Float64(file_dist_km), Float64(dist_km); rtol=0, atol=eps(Float64) * max(abs(dist_km), 1)) ||
            throw(ArgumentError("particle file was generated with dist_km=$file_dist_km, not $dist_km"))
    end
    if haskey(ds.attrib, _PARTICLE_ORDER_ATTRIBUTE)
        String(ds.attrib[_PARTICLE_ORDER_ATTRIBUTE]) == _PARTICLE_ORDER ||
            throw(ArgumentError("particle file has unsupported FTLE particle order"))
    end
    return nothing
end

function _initial_position_tolerance(expected_plonds, expected_platds)
    max_offset = 0.0
    for i in firstindex(expected_plonds):4:lastindex(expected_plonds)
        max_offset = max(
            max_offset,
            abs(_wrapped_lon_diff(expected_plonds[i], expected_plonds[i + 1])) / 2,
            abs(expected_platds[i + 2] - expected_platds[i + 3]) / 2,
        )
    end
    return max(1.2e-2, 0.05 * max_offset)
end

function _validate_initial_positions(ds, grid_or_spectral_grid, dist_km)
    grid_or_spectral_grid isa Integer &&
        throw(ArgumentError("initial position validation needs a grid; pass validate_initial_positions=false for npoint-only files"))

    expected_plonds, expected_platds = initial_FTLE_particle_positions(grid_or_spectral_grid, dist_km)
    actual_plonds = ds["lon"][:, 1]
    actual_platds = ds["lat"][:, 1]
    tolerance = _initial_position_tolerance(expected_plonds, expected_platds)

    for i in eachindex(expected_plonds)
        lon_error = abs(_wrapped_lon_diff(Float64(actual_plonds[i]), expected_plonds[i]))
        lat_error = abs(Float64(actual_platds[i]) - expected_platds[i])
        max(lon_error, lat_error) <= tolerance ||
            throw(ArgumentError("particle file initial positions do not match the FTLE stencil"))
    end

    return nothing
end

function _validate_particle_file(ds, grid_or_spectral_grid, dist_km; validate_initial_positions)
    haskey(ds.dim, "particle") || throw(ArgumentError("particle file must contain a particle dimension"))
    haskey(ds.dim, "time") || throw(ArgumentError("particle file must contain a time dimension"))
    for name in ("time", "lon", "lat")
        haskey(ds, name) || throw(ArgumentError("particle file must contain variable $name"))
    end

    npoints = _grid_npoints(grid_or_spectral_grid)
    nparticles = ds.dim["particle"]
    ntimes = ds.dim["time"]
    expected_particles = 4 * npoints
    nparticles == expected_particles ||
        throw(DimensionMismatch("particle file has $nparticles particles, expected $expected_particles"))
    ntimes > 0 || throw(ArgumentError("particle file time dimension is empty"))

    dimnames(ds["time"]) == ("time",) ||
        throw(DimensionMismatch("particle file time variable must have dimension (time)"))
    for name in ("lon", "lat")
        dimnames(ds[name]) == ("particle", "time") ||
            throw(DimensionMismatch("particle file $name variable must have dimensions (particle, time)"))
        Tuple(dimsize(ds[name])) == (nparticles, ntimes) ||
            throw(DimensionMismatch("particle file $name variable has the wrong size"))
    end

    _validate_particle_file_metadata(ds, dist_km)
    validate_initial_positions && _validate_initial_positions(ds, grid_or_spectral_grid, dist_km)
    return nothing
end

function _write_particle_file_metadata(path; dist_km)
    ds = NCDataset(path, "a")
    try
        ds.attrib[_DIST_KM_ATTRIBUTE] = Float64(dist_km)
        ds.attrib[_PARTICLE_ORDER_ATTRIBUTE] = _PARTICLE_ORDER
    finally
        close(ds)
    end
    return path
end

"""
    FTLE_from_particle_file!(ftle, B, path, grid_or_spectral_grid, dist_km; time_indices = :, validate_initial_positions = true)

Post-process a SpeedyWeather `ParticleTracker` NetCDF file into caller-provided
output and work arrays. Returns `ftle, selected_time_hours`.
"""
function FTLE_from_particle_file!(
    ftle,
    B,
    path::AbstractString,
    grid_or_spectral_grid,
    dist_km;
    time_indices=Colon(),
    validate_initial_positions=true,
)
    ds = NCDataset(path, "r")
    try
        _validate_particle_file(ds, grid_or_spectral_grid, dist_km; validate_initial_positions)
        time_hours = _particle_file_time_hours(ds)
        selected = _checked_time_indices(time_indices, time_hours)
        FTLE_from_particles!(ftle, B, ds["lon"], ds["lat"], time_hours, grid_or_spectral_grid, dist_km; time_indices=selected)
        return ftle, _selected_time_hours(time_hours, selected)
    finally
        close(ds)
    end
end

"""
    FTLE_from_particle_file(path, grid_or_spectral_grid, dist_km; time_indices = :, validate_initial_positions = true)

Allocate and compute FTLE from a SpeedyWeather `ParticleTracker` NetCDF file.
Returns `ftle, selected_time_hours`.
"""
function FTLE_from_particle_file(path::AbstractString, grid_or_spectral_grid, dist_km; time_indices=Colon(), validate_initial_positions=true)
    ds = NCDataset(path, "r")
    try
        _validate_particle_file(ds, grid_or_spectral_grid, dist_km; validate_initial_positions)
        time_hours = _particle_file_time_hours(ds)
        selected = _checked_time_indices(time_indices, time_hours)
        npoints = _grid_npoints(grid_or_spectral_grid)
        ftle = Array{Float64}(undef, npoints, length(selected))
        B = Array{Float64}(undef, 2, 2, npoints)
        FTLE_from_particles!(ftle, B, ds["lon"], ds["lat"], time_hours, npoints, dist_km; time_indices=selected)
        return ftle, _selected_time_hours(time_hours, selected)
    finally
        close(ds)
    end
end

"""
    FTLEResult(ftle, spectral_grid, time_hours; dist_km, backwards, ...)

Small named container returned by `get_FTLE(...; return_result = true)`.
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
    backwards=false,
    direction=backwards ? :negative : :positive,
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
        particle_file_path,
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

"""
    final_ftle(result::FTLEResult)

Return the final selected FTLE column from `result`.
"""
final_ftle(result::FTLEResult) = result.ftle[:, end]

function _nearest_time_index(times, time_hour)
    isempty(times) && throw(ArgumentError("no time_hours are available"))
    isfinite(time_hour) || throw(ArgumentError("time_hour must be finite"))
    distances = abs.(times .- time_hour)
    return first(eachindex(times)) + argmin(distances) - 1
end

function _resolve_time_index(ftle::AbstractMatrix; time_index, time_hour, time_hours)
    time_index !== nothing && time_hour !== nothing &&
        throw(ArgumentError("pass either time_index or time_hour, not both"))
    if time_hour !== nothing
        time_hours === nothing && throw(ArgumentError("time_hour requires time_hours"))
        length(time_hours) == size(ftle, 2) ||
            throw(DimensionMismatch("time_hours must contain one value per FTLE column"))
        time_index = _nearest_time_index(time_hours, time_hour)
    elseif time_index === nothing
        time_index = size(ftle, 2)
    end
    1 <= time_index <= size(ftle, 2) || throw(BoundsError(ftle, (:, time_index)))
    return time_index
end

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
    result.spectral_grid === nothing &&
        throw(ArgumentError("ftle_field requires result.spectral_grid"))
    if time_hour !== nothing
        isequal(time_indices, Colon()) ||
            throw(ArgumentError("pass either time_indices or time_hour, not both"))
        index = _resolve_time_index(result.ftle; time_index=nothing, time_hour, time_hours=result.time_hours)
        return ftle_field(view(result.ftle, :, index), result.spectral_grid)
    end

    selected = _checked_time_indices(time_indices, result.time_hours)
    length(selected) == 1 && return ftle_field(view(result.ftle, :, first(selected)), result.spectral_grid)
    return ftle_field(view(result.ftle, :, selected), result.spectral_grid)
end

"""
    final_ftle_field(result::FTLEResult)

Return the final selected FTLE column as a `RingGrids.Field`.
"""
final_ftle_field(result::FTLEResult) = ftle_field(final_ftle(result), result.spectral_grid)

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

function Base.show(io::IO, result::FTLEResult)
    print(io, "FTLEResult(", size(result.ftle, 1), " points, ", size(result.ftle, 2), " times, direction=", result.direction, ")")
end

const _SPEEDYWEATHER_SOURCE_INSTALL_MESSAGE =
    "SpeedyWeatherFTLE requires the mk/lyapunov2 SpeedyWeather sources because registered releases do not yet expose the needed particle-advection API."

function _particle_advection_2d(spectral_grid; nparticles, backwards, every_n_time_steps)
    try
        return ParticleAdvection2D(spectral_grid; nparticles, backwards, every_n_time_steps)
    catch err
        err isa MethodError && throw(ArgumentError(_SPEEDYWEATHER_SOURCE_INSTALL_MESSAGE))
        rethrow()
    end
end

_period_hours(hours) = Second(round(Int, 3600 * _check_positive(hours, "hours")))
_period_days(days) = Second(round(Int, 86_400 * _check_positive(days, "simulation_days")))

function _spectral_grid_for(field::Field)
    grid = field.grid
    grid_type = typeof(grid)
    rings = getproperty(grid, :rings)
    guesses = unique(Int[floor(Int, 2 * length(rings) / 3), ceil(Int, 2 * length(rings) / 3)])
    for trunc in guesses
        spectral_grid = SpectralGrid(nlayers=1, trunc=trunc, Grid=grid_type)
        spectral_grid.grid == grid && return spectral_grid
    end
    throw(ArgumentError("could not infer a SpectralGrid for the supplied field grid"))
end

"""
    get_FTLE(u::Field, v::Field; kwargs...)

Run frozen-flow SpeedyWeather particle tracking for zonal and meridional
velocity fields, then return `(ftle, spectral_grid, time_hours)`.
"""
function get_FTLE(
    u::Field,
    v::Field;
    simulation_days=10,
    dist_km=10,
    backwards=false,
    dynamics=false,
    rint_hours=3,
    particle_advection_every_n_time_steps=6,
    particle_tracker_keepbits=15,
    particle_tracker_path="",
    particle_tracker_filename="particles.nc",
    keep_particle_file=false,
    return_particle_file_path=false,
    return_result=false,
    time_indices=Colon(),
)
    u.grid == v.grid || throw(ArgumentError("u and v must use the same grid"))
    dist_km = _check_dist_km(dist_km)
    dynamics && throw(ArgumentError("get_FTLE supports frozen prescribed fields only; pass dynamics=false"))
    particle_advection_every_n_time_steps >= 1 ||
        throw(ArgumentError("particle_advection_every_n_time_steps must be at least 1"))
    rint_hours = _check_positive(rint_hours, "rint_hours")
    particle_tracker_keepbits >= 1 ||
        throw(ArgumentError("particle_tracker_keepbits must be positive"))

    spectral_grid = _spectral_grid_for(u)
    nparticles = 4 * spectral_grid.npoints
    particle_advection = _particle_advection_2d(
        spectral_grid;
        nparticles,
        backwards,
        every_n_time_steps=particle_advection_every_n_time_steps,
    )

    model = BarotropicModel(spectral_grid; dynamics=false, particle_advection)
    simulation = initialize!(model)
    model.callbacks[:particle_tracker] = ParticleTracker(
        spectral_grid;
        schedule=Schedule(every=_period_hours(rint_hours)),
        keepbits=particle_tracker_keepbits,
        path=particle_tracker_path,
        filename=particle_tracker_filename,
    )

    londs, latds = RingGrids.get_londlatds(u.grid)
    perturb_positions_FTLE(simulation.variables.prognostic.particles, londs, latds, dist_km)

    SpeedyWeather.initialize!(simulation; period=_period_days(simulation_days))
    simulation.variables.grid.u[:, 1, 1] .= u
    simulation.variables.grid.v[:, 1, 1] .= v
    SpeedyWeather.initialize!(simulation.variables, simulation.variables.prognostic.particles, model)
    SpeedyWeather.time_stepping!(simulation)
    SpeedyWeather.finalize!(simulation)

    tracker = model.callbacks[:particle_tracker]
    path = joinpath(tracker.path == "" ? model.output.run_path : tracker.path, tracker.filename)
    keep_file = keep_particle_file || return_particle_file_path

    try
        _write_particle_file_metadata(path; dist_km)
        ftle, time_hours = FTLE_from_particle_file(path, spectral_grid, dist_km; time_indices)
        if return_result
            return FTLEResult(
                ftle,
                spectral_grid,
                time_hours;
                particle_file_path=keep_file ? path : nothing,
                dist_km,
                backwards,
                dynamics=false,
                rint_hours,
            )
        elseif return_particle_file_path
            return ftle, spectral_grid, time_hours, path
        else
            return ftle, spectral_grid, time_hours
        end
    finally
        keep_file || rm(path; force=true)
    end
end

function _reject_backwards(kwargs, name)
    :backwards in keys(kwargs) &&
        throw(ArgumentError("$name fixes the FTLE time direction; call get_FTLE to choose it explicitly"))
    return nothing
end

"""
    positive_FTLE(u::Field, v::Field; kwargs...)

Compute forward-time FTLE for a frozen prescribed velocity field by calling
[`get_FTLE`](@ref) with `backwards = false`.
"""
function positive_FTLE(u::Field, v::Field; kwargs...)
    _reject_backwards(kwargs, "positive_FTLE")
    return get_FTLE(u, v; backwards=false, kwargs...)
end

"""
    negative_FTLE(u::Field, v::Field; kwargs...)

Compute backward-time FTLE for a frozen prescribed velocity field by calling
[`get_FTLE`](@ref) with `backwards = true`.
"""
function negative_FTLE(u::Field, v::Field; kwargs...)
    _reject_backwards(kwargs, "negative_FTLE")
    return get_FTLE(u, v; backwards=true, kwargs...)
end

end
