@inline _wrapped_lon_diff(lond1, lond2) = mod(lond1 - lond2 + 180, 360) - 180

const _FTLE_PARTICLE_FILE_METADATA_PREFIX = "SpeedyWeatherFTLE_"
const _FTLE_PARTICLE_FILE_DIST_KM_ATTRIBUTE = _FTLE_PARTICLE_FILE_METADATA_PREFIX * "dist_km"
const _FTLE_PARTICLE_FILE_ORDER_ATTRIBUTE = _FTLE_PARTICLE_FILE_METADATA_PREFIX * "particle_order"
const _FTLE_PARTICLE_FILE_KEEPBITS_ATTRIBUTE = _FTLE_PARTICLE_FILE_METADATA_PREFIX * "particle_tracker_keepbits"
const _FTLE_PARTICLE_ORDER = "east,west,north,south"

# ParticleTracker stores lon/lat as Float32 and may quantize mantissa bits.
# This tolerance covers normal Float32/keepbits=15 output while rejecting common
# 10 km vs 20 km stencil mistakes. Coarsely quantized legacy files should carry
# metadata or opt out after external validation.
const _INITIAL_POSITION_MIN_ATOL_DEGREES = 1.2e-2
const _INITIAL_POSITION_RTOL = 5.0e-2

function _check_dist_km(dist_km)
    dist_km isa Real && isfinite(dist_km) && dist_km > 0 ||
        throw(ArgumentError("dist_km must be finite and positive"))
    return dist_km
end

function _check_positive_hours(value, name)
    value isa Real && isfinite(value) && value > 0 ||
        throw(ArgumentError("$name must be finite and positive"))
    return value
end

function _duration_magnitude(time_value, name)
    time_value isa Real && isfinite(time_value) ||
        throw(ArgumentError("$name must be finite"))
    return abs(float(time_value))
end

function _check_time_hours(time_hours)
    for i in eachindex(time_hours)
        time_hours[i] isa Real && isfinite(time_hours[i]) ||
            throw(ArgumentError("time_hours[$i] must be real and finite"))
    end
    return nothing
end

function displacement_gradient_matrix_central!(B, plonds, platds, dist_km)
    """
    Compute the flow-map Jacobian given particle positions at a fixed time
    Uses a central difference scheme to do this, with four points per grid cell

    Inputs:
        plonds: longitudes of particles
        platds: latitudes of particles
        dist_km: perturbation in position applied before starting simulation

    Outputs:
        B: (2,2,N) array where N is the number of grid points in the simulation.
        B[:,:,k] is the flow-map Jacobian for the k'th grid point.
    """

    _check_dist_km(dist_km)
    length(plonds) == length(platds) || throw(DimensionMismatch("plonds and platds must have the same length"))
    length(plonds) % 4 == 0 || throw(ArgumentError("particle position vectors must contain four particles per grid point"))

    Ngpoints = length(plonds) ÷ 4 # Number of grid points
    size(B) == (2, 2, Ngpoints) || throw(DimensionMismatch("B must have size (2, 2, $Ngpoints)"))

    dfac = Re / (dist_km * 1000) / 2

    @inbounds for i in 1:Ngpoints
        p = 4i

        # Derivative of x w.r.t. X
        B[1, 1, i] = deg2rad(_wrapped_lon_diff(plonds[p - 3], plonds[p - 2])) *
                     cosd((platds[p - 3] + platds[p - 2]) / 2) * dfac
        # Derivative of y w.r.t. X
        B[2, 1, i] = deg2rad(platds[p - 3] - platds[p - 2]) * dfac
        # Derivative of x w.r.t. Y
        B[1, 2, i] = deg2rad(_wrapped_lon_diff(plonds[p - 1], plonds[p])) *
                     cosd((platds[p - 1] + platds[p]) / 2) * dfac
        # Derivative of y w.r.t. Y
        B[2, 2, i] = deg2rad(platds[p - 1] - platds[p]) * dfac
    end

    return B
end

function displacement_gradient_matrix_central(plonds, platds, dist_km)
    Ngpoints = length(plonds) ÷ 4 # Number of grid points
    B = Array{Float64}(undef, 2, 2, Ngpoints)
    displacement_gradient_matrix_central!(B, plonds, platds, dist_km)
    return B
end

@inline function _largest_cauchy_green_eigenvalue(a, b, c, d)
    C11 = a*a + c*c
    C12 = a*b + c*d
    C22 = b*b + d*d
    return (C11 + C22 + sqrt((C11 - C22)^2 + 4C12^2)) / 2
end

function FTLE_over_grid!(FTLE_grid, B, T)
    """
    Compute FTLE over a grid

    Inputs:
        B: flow-map Jacobian/deformation gradient
        T: elapsed integration duration. Signed values are accepted, and the
           FTLE rate is computed using `abs(T)`.

    Outputs: 
        FTLE_grid: FTLE at each grid point. Units are 1 / [T].
        If T is zero, the FTLE is undefined and NaN is written. Non-finite
        durations throw `ArgumentError`.
    """

    Ngpoints = size(B, 3) # Number of grid points
    length(FTLE_grid) == Ngpoints || throw(DimensionMismatch("FTLE_grid must have length $Ngpoints"))

    duration = _duration_magnitude(T, "T")
    if iszero(duration)
        fill!(FTLE_grid, NaN)
        return FTLE_grid
    end

    twoT = 2 * duration
    for k in 1:Ngpoints
        # Largest eigenvalue of the 2x2 right Cauchy-Green tensor B'B.
        @inbounds lmax = _largest_cauchy_green_eigenvalue(B[1, 1, k], B[1, 2, k], B[2, 1, k], B[2, 2, k])
        # FTLE - in units of 1 / [T]
        @inbounds FTLE_grid[k] = log(lmax) / twoT
    end

    return FTLE_grid
end

function FTLE_over_grid(B, T)
    FTLE_grid = Vector{Float64}(undef, size(B, 3))
    FTLE_over_grid!(FTLE_grid, B, T)
    return FTLE_grid
end

_particle_column(A::Union{NCDatasets.Variable, NCDatasets.CFVariable}, tindex) = A[:, tindex]
_particle_column(A::AbstractArray, tindex) = view(A, :, tindex)
_particle_column(A, tindex) = A[:, tindex]

function _check_time_index(tindex::Integer, n_times)
    tindex isa Bool &&
        throw(ArgumentError("time_indices must be integer indices, not Bool values"))
    1 <= tindex <= n_times || throw(BoundsError(1:n_times, tindex))
    return nothing
end

_checked_time_indices(time_indices, time_hours::AbstractVector{<:Real}) =
    _checked_time_indices(time_indices, length(time_hours), time_hours)

_checked_time_indices(::Colon, n_times, time_hours) = Base.OneTo(n_times)

function _checked_time_indices(time_indices::Symbol, n_times, time_hours)
    if time_indices === :all
        return Base.OneTo(n_times)
    elseif time_indices === :first
        return (1,)
    elseif time_indices in (:last, :final)
        return (n_times,)
    elseif time_indices === :nonzero
        return findall(t -> isfinite(t) && !iszero(t), time_hours)
    elseif time_indices === :positive
        return findall(t -> isfinite(t) && t > 0, time_hours)
    else
        throw(ArgumentError("unsupported time_indices selector :$time_indices; use :, :all, :first, :last, :final, :nonzero, :positive, an integer index, or integer indices"))
    end
end

function _checked_time_indices(tindex::Integer, n_times, time_hours)
    _check_time_index(tindex, n_times)
    return (tindex,)
end

function _checked_time_indices(time_indices::Real, n_times, time_hours)
    throw(ArgumentError("time_indices must be :, a supported Symbol selector, an integer index, or integer indices"))
end

function _checked_time_indices(time_indices::AbstractVector{Bool}, n_times, time_hours)
    length(time_indices) == n_times ||
        throw(DimensionMismatch("boolean time_indices mask has length $(length(time_indices)), but time_hours has length $n_times"))
    return findall(time_indices)
end

function _checked_time_indices(time_indices, n_times, time_hours)
    for tindex in time_indices
        tindex isa Integer ||
            throw(ArgumentError("time_indices must contain only integer indices"))
        _check_time_index(tindex, n_times)
    end
    return time_indices
end

_selected_time_hours(time_hours, time_indices) = [Float64(time_hours[tindex]) for tindex in time_indices]

function _check_particle_position_column(plonds, platds, tindex)
    length(plonds) == length(platds) ||
        throw(DimensionMismatch("particle longitude and latitude columns at time index $tindex must have the same length"))

    for i in eachindex(plonds, platds)
        lon = plonds[i]
        lat = platds[i]
        (ismissing(lon) || ismissing(lat)) &&
            throw(ArgumentError("particle positions at time index $tindex contain missing values"))
        lon isa Real && lat isa Real ||
            throw(ArgumentError("particle positions at time index $tindex must be real-valued"))
        isfinite(lon) && isfinite(lat) ||
            throw(ArgumentError("particle positions at time index $tindex contain non-finite values"))
    end

    return nothing
end

"""
    FTLE_from_particles!(
        FTLE_grid_time,
        B,
        plonds_time,
        platds_time,
        time_hours,
        grid_or_npoints,
        dist_km;
        time_indices = :
    )

Compute FTLE from in-memory particle trajectories, writing into caller-provided
arrays.

`plonds_time` and `platds_time` must have dimensions `(particle, time)`, with
four particles per FTLE grid point in east, west, north, south order. `time_hours`
contains the corresponding sampling times in hours. `grid_or_npoints` may be a
spatial grid, a `SpectralGrid`, or the number of FTLE grid points. `dist_km` is
the initial perturbation distance used to release the particles.

`FTLE_grid_time` must have dimensions `(grid point, selected time)`, where the
selected columns are determined by `time_indices`. `B` is a reusable
deformation-gradient work array with dimensions `(2, 2, grid point)`.

Supported `time_indices` values are `:`, `:all`, `:first`, `:last`, `:final`,
`:nonzero`, `:positive`, an integer index, integer-index iterables, or a boolean
mask. Signed sampling times are treated as elapsed durations with magnitude
`abs(time_hours)`. Zero-duration samples are filled with `NaN` because FTLE is
undefined at `t = 0`; non-finite sampling times throw `ArgumentError`.

Returns `FTLE_grid_time`.
"""
function FTLE_from_particles!(
    FTLE_grid_time,
    B,
    plonds_time,
    platds_time,
    time_hours::AbstractVector{<:Real},
    grid_or_npoints,
    dist_km;
    time_indices=Colon(),
)
    """
    Compute FTLE from particle trajectories, reusing caller-provided output and work arrays.

    Inputs:
        FTLE_grid_time: output matrix with dimensions (grid point, time)
        B: deformation-gradient work array with dimensions (2, 2, grid point)
        plonds_time: longitude matrix with dimensions (particle, time)
        platds_time: latitude matrix with dimensions (particle, time)
        time_hours: sampling times in hours
        grid_or_npoints: spatial grid, SpectralGrid, or number of FTLE grid points
        dist_km: initial FTLE particle perturbation in km
        time_indices: optional time columns to process. Supports `:`, `:all`,
            `:first`, `:last`/`:final`, `:nonzero`/`:positive`, integer
            indices, integer-index iterables, and boolean masks.

    Output:
        FTLE_grid_time, modified in place.
    """

    _check_dist_km(dist_km)
    size(plonds_time) == size(platds_time) ||
        throw(DimensionMismatch("plonds_time and platds_time must have the same size"))

    npoints = _grid_npoints(grid_or_npoints)
    npoints > 0 || throw(ArgumentError("number of grid points must be positive"))

    expected_nparticles = 4 * npoints
    actual_nparticles = size(plonds_time, 1)
    actual_nparticles == expected_nparticles ||
        throw(DimensionMismatch("particle trajectories have $actual_nparticles particles, expected $expected_nparticles"))

    n_times = length(time_hours)
    _check_time_hours(time_hours)
    size(plonds_time, 2) == n_times ||
        throw(DimensionMismatch("particle trajectories have $(size(plonds_time, 2)) time samples, but time_hours has length $n_times"))
    selected_time_indices = _checked_time_indices(time_indices, time_hours)
    n_selected_times = length(selected_time_indices)
    size(FTLE_grid_time) == (npoints, n_selected_times) ||
        throw(DimensionMismatch("FTLE_grid_time must have size ($npoints, $n_selected_times)"))
    size(B) == (2, 2, npoints) ||
        throw(DimensionMismatch("B must have size (2, 2, $npoints)"))

    for (out_index, tindex) in enumerate(selected_time_indices)
        thour = time_hours[tindex]
        plonds = _particle_column(plonds_time, tindex)
        platds = _particle_column(platds_time, tindex)

        _check_particle_position_column(plonds, platds, tindex)
        displacement_gradient_matrix_central!(B, plonds, platds, dist_km)
        FTLE_over_grid!(view(FTLE_grid_time, :, out_index), B, thour)
    end

    return FTLE_grid_time
end

"""
    FTLE_from_particles(
        plonds_time,
        platds_time,
        time_hours,
        grid_or_npoints,
        dist_km;
        time_indices = :
    )

Compute FTLE from in-memory particle trajectories and allocate the output
arrays.

This is the allocating companion to [`FTLE_from_particles!`](@ref). It returns
`FTLE_grid_time, selected_time_hours`, where `FTLE_grid_time` has dimensions
`(grid point, selected time)`.
"""
function FTLE_from_particles(
    plonds_time,
    platds_time,
    time_hours::AbstractVector{<:Real},
    grid_or_npoints,
    dist_km;
    time_indices=Colon(),
)
    """
    Compute FTLE from in-memory particle trajectories.

    Inputs:
        plonds_time: longitude matrix with dimensions (particle, time)
        platds_time: latitude matrix with dimensions (particle, time)
        time_hours: sampling times in hours
        grid_or_npoints: spatial grid, SpectralGrid, or number of FTLE grid points
        dist_km: initial FTLE particle perturbation in km
        time_indices: optional time columns to process. Supports `:`, `:all`,
            `:first`, `:last`/`:final`, `:nonzero`/`:positive`, integer
            indices, integer-index iterables, and boolean masks.

    Outputs:
        FTLE_grid_time: FTLE matrix with one row per grid point and one column per recorded time
        time_hours: selected sampling times in hours, converted to Float64

    Signed sampling times are treated as elapsed durations with magnitude
    `abs(time_hours)`. Zero-duration samples are returned as NaN because FTLE is
    undefined at T = 0. Non-finite sampling times throw `ArgumentError`.
    """

    npoints = _grid_npoints(grid_or_npoints)
    n_times = length(time_hours)
    _check_time_hours(time_hours)
    selected_time_indices = _checked_time_indices(time_indices, time_hours)
    time_hours_float = _selected_time_hours(time_hours, selected_time_indices)
    FTLE_grid_time = Array{Float64}(undef, npoints, length(selected_time_indices))
    B = Array{Float64}(undef, 2, 2, npoints)

    FTLE_from_particles!(
        FTLE_grid_time,
        B,
        plonds_time,
        platds_time,
        time_hours,
        npoints,
        dist_km;
        time_indices=selected_time_indices,
    )
    return FTLE_grid_time, time_hours_float
end

function _particle_file_time_hours(particles_ds)
    time_vec = particles_ds["time"][:]
    isempty(time_vec) &&
        throw(ArgumentError("particle file time variable must contain at least one sample"))
    first_time = time_vec[begin]
    if first_time isa Real
        _check_time_hours(time_vec)
        return Float64.(time_vec .- first_time)
    else
        try
            return Float64.((time_vec .- first_time) ./ Hour(1)) # Hours since release time
        catch
            throw(ArgumentError("particle file time variable must contain datetime samples or real-valued elapsed times"))
        end
    end
end

function _expected_initial_particle_positions(grid_or_spectral_grid, dist_km)
    return initial_FTLE_particle_positions(grid_or_spectral_grid, dist_km)
end

function _expected_initial_particle_positions(npoints::Integer, dist_km)
    throw(ArgumentError("initial-position validation requires a grid or SpectralGrid; pass validate_initial_positions=false to skip this compatibility check"))
end

function _particle_file_attribute(particles_ds, name)
    haskey(particles_ds.attrib, name) || return nothing
    return particles_ds.attrib[name]
end

function _write_FTLE_particle_file_metadata(path; dist_km, particle_tracker_keepbits)
    particles_ds = NCDataset(path, "a")
    try
        particles_ds.attrib[_FTLE_PARTICLE_FILE_DIST_KM_ATTRIBUTE] = Float64(dist_km)
        particles_ds.attrib[_FTLE_PARTICLE_FILE_ORDER_ATTRIBUTE] = _FTLE_PARTICLE_ORDER
        particles_ds.attrib[_FTLE_PARTICLE_FILE_KEEPBITS_ATTRIBUTE] = Int(particle_tracker_keepbits)
    finally
        close(particles_ds)
    end
    return path
end

function _validate_particle_file_metadata(particles_ds, dist_km)
    metadata_dist_km = _particle_file_attribute(particles_ds, _FTLE_PARTICLE_FILE_DIST_KM_ATTRIBUTE)
    if metadata_dist_km !== nothing
        metadata_dist_km isa Real ||
            throw(ArgumentError("particle file FTLE metadata dist_km must be real-valued"))
        isapprox(Float64(metadata_dist_km), Float64(dist_km); rtol=0, atol=eps(Float64) * max(abs(Float64(dist_km)), 1.0)) ||
            throw(ArgumentError("particle file was generated with dist_km=$(metadata_dist_km), but post-processing requested dist_km=$(dist_km)"))
    end

    metadata_order = _particle_file_attribute(particles_ds, _FTLE_PARTICLE_FILE_ORDER_ATTRIBUTE)
    if metadata_order !== nothing && String(metadata_order) != _FTLE_PARTICLE_ORDER
        throw(ArgumentError("particle file FTLE metadata has particle_order=$(metadata_order), expected $(_FTLE_PARTICLE_ORDER)"))
    end

    return nothing
end

function _validate_particle_file_structure(particles_ds, grid_or_spectral_grid)
    haskey(particles_ds.dim, "particle") ||
        throw(ArgumentError("particle file must contain a particle dimension"))
    haskey(particles_ds.dim, "time") ||
        throw(ArgumentError("particle file must contain a time dimension"))

    npoints = _grid_npoints(grid_or_spectral_grid)
    expected_nparticles = 4 * npoints
    actual_nparticles = particles_ds.dim["particle"]
    actual_nparticles == expected_nparticles ||
        throw(DimensionMismatch("particle file has $actual_nparticles particles, expected $expected_nparticles"))

    n_times = particles_ds.dim["time"]
    n_times > 0 ||
        throw(ArgumentError("particle file time dimension must contain at least one sample"))

    for variable_name in ("time", "lon", "lat")
        haskey(particles_ds, variable_name) ||
            throw(ArgumentError("particle file must contain a $variable_name variable"))
    end

    time_dims = dimnames(particles_ds["time"])
    time_dims == ("time",) ||
        throw(DimensionMismatch("particle file time variable must have dimensions (time), got $time_dims"))
    length(particles_ds["time"]) == n_times ||
        throw(DimensionMismatch("particle file time variable has length $(length(particles_ds["time"])), expected $n_times"))

    expected_position_dims = ("particle", "time")
    for variable_name in ("lon", "lat")
        variable = particles_ds[variable_name]
        variable_dims = dimnames(variable)
        variable_dims == expected_position_dims ||
            throw(DimensionMismatch("particle file $variable_name variable must have dimensions $expected_position_dims, got $variable_dims"))
        variable_size = Tuple(dimsize(variable))
        expected_size = (actual_nparticles, n_times)
        variable_size == expected_size ||
            throw(DimensionMismatch("particle file $variable_name variable has size $variable_size, expected $expected_size"))
    end

    return npoints, n_times
end

function _initial_position_atol_degrees(expected_plonds, expected_platds)
    max_offset = 0.0
    for i in firstindex(expected_plonds):4:lastindex(expected_plonds)
        max_offset = max(
            max_offset,
            abs(_wrapped_lon_diff(expected_plonds[i], expected_plonds[i + 1])) / 2,
            abs(expected_platds[i + 2] - expected_platds[i + 3]) / 2,
        )
    end
    return max(_INITIAL_POSITION_MIN_ATOL_DEGREES, _INITIAL_POSITION_RTOL * max_offset)
end

function _validate_particle_file_initial_positions(particles_ds, grid_or_spectral_grid, dist_km)
    expected_plonds, expected_platds = _expected_initial_particle_positions(grid_or_spectral_grid, dist_km)
    actual_plonds = particles_ds["lon"][:, 1]
    actual_platds = particles_ds["lat"][:, 1]
    length(actual_plonds) == length(expected_plonds) ||
        throw(DimensionMismatch("particle file initial longitude column has length $(length(actual_plonds)), expected $(length(expected_plonds))"))

    max_lon_error = 0.0
    max_lat_error = 0.0
    for i in eachindex(expected_plonds)
        actual_lon = actual_plonds[i]
        actual_lat = actual_platds[i]
        (ismissing(actual_lon) || ismissing(actual_lat)) &&
            throw(ArgumentError("particle file initial positions contain missing values"))

        lon_error = abs(_wrapped_lon_diff(Float64(actual_lon), expected_plonds[i]))
        lat_error = abs(Float64(actual_lat) - expected_platds[i])
        (isfinite(lon_error) && isfinite(lat_error)) ||
            throw(ArgumentError("particle file initial positions contain non-finite values"))
        max_lon_error = max(max_lon_error, lon_error)
        max_lat_error = max(max_lat_error, lat_error)
    end

    atol_degrees = _initial_position_atol_degrees(expected_plonds, expected_platds)
    if max_lon_error > atol_degrees || max_lat_error > atol_degrees
        throw(ArgumentError("particle file initial positions do not match the FTLE east/west/north/south stencil for this grid and dist_km (max longitude error $(max_lon_error)°, max latitude error $(max_lat_error)°, tolerance $(atol_degrees)°); if this is a coarsely quantized legacy file, regenerate it with SpeedyWeatherFTLE metadata or pass validate_initial_positions=false only after external validation"))
    end

    return nothing
end

function _validate_particle_file(particles_ds, grid_or_spectral_grid, dist_km; validate_initial_positions=true)
    _check_dist_km(dist_km)
    _validate_particle_file_structure(particles_ds, grid_or_spectral_grid)
    _validate_particle_file_metadata(particles_ds, dist_km)

    if validate_initial_positions
        _validate_particle_file_initial_positions(particles_ds, grid_or_spectral_grid, dist_km)
    end

    return nothing
end

"""
    FTLE_from_particle_file!(
        FTLE_grid_time,
        B,
        path,
        grid_or_spectral_grid,
        dist_km;
        time_indices = :,
        validate_initial_positions = true
    )

Compute FTLE from a SpeedyWeather `ParticleTracker` NetCDF file, writing into
caller-provided arrays.

`path` must point to a particle-tracker file with longitude and latitude
variables named `lon` and `lat`. By default, the file must contain the
canonical FTLE-compatible initial stencil: four particles per FTLE grid point in
east, west, north, south order, matching `grid_or_spectral_grid` and `dist_km`.
Pass `validate_initial_positions = false` only when intentionally processing a
pre-validated compatible file whose initial positions cannot be checked.
`grid_or_spectral_grid` may be the spatial grid or the `SpectralGrid` used for
the tracking run.

The output and work arrays have the same requirements as
[`FTLE_from_particles!`](@ref). Returns `FTLE_grid_time, selected_time_hours`.
"""
function FTLE_from_particle_file!(
    FTLE_grid_time,
    B,
    path::AbstractString,
    grid_or_spectral_grid,
    dist_km;
    time_indices=Colon(),
    validate_initial_positions=true,
)
    """
    Compute FTLE from a SpeedyWeather `ParticleTracker` NetCDF file, reusing
    caller-provided output and work arrays.

    Inputs:
        FTLE_grid_time: output matrix with dimensions (grid point, time)
        B: deformation-gradient work array with dimensions (2, 2, grid point)
        path: path to the particle-tracker NetCDF file
        grid_or_spectral_grid: spatial grid, or the SpectralGrid used for tracking
        dist_km: initial FTLE particle perturbation in km
        validate_initial_positions: validate the first particle positions
            against the canonical FTLE stencil before post-processing
        time_indices: optional time columns to process. Supports `:`, `:all`,
            `:first`, `:last`/`:final`, `:nonzero`/`:positive`, integer
            indices, integer-index iterables, and boolean masks.

    Outputs:
        FTLE_grid_time: modified in place
        time_hours: selected sampling times in hours

    Signed sampling times are treated as elapsed durations with magnitude
    `abs(time_hours)`. Zero-duration samples are returned as NaN because FTLE is
    undefined at T = 0. Non-finite sampling times throw `ArgumentError`.
    """

    npoints = _grid_npoints(grid_or_spectral_grid)
    particles_ds = NCDataset(path, "r")

    try
        _validate_particle_file(
            particles_ds,
            grid_or_spectral_grid,
            dist_km;
            validate_initial_positions,
        )
        time_hours = _particle_file_time_hours(particles_ds)
        _check_time_hours(time_hours)
        selected_time_indices = _checked_time_indices(time_indices, time_hours)
        FTLE_from_particles!(
            FTLE_grid_time,
            B,
            particles_ds["lon"],
            particles_ds["lat"],
            time_hours,
            npoints,
            dist_km;
            time_indices=selected_time_indices,
        )
        return FTLE_grid_time, _selected_time_hours(time_hours, selected_time_indices)
    finally
        close(particles_ds)
    end
end

"""
    FTLE_from_particle_file(
        path,
        grid_or_spectral_grid,
        dist_km;
        time_indices = :,
        validate_initial_positions = true
    )

Compute FTLE from a SpeedyWeather `ParticleTracker` NetCDF file and allocate the
output arrays.

This is the allocating companion to [`FTLE_from_particle_file!`](@ref). It
returns `FTLE_grid_time, selected_time_hours`, where `FTLE_grid_time` has
dimensions `(grid point, selected time)`.
"""
function FTLE_from_particle_file(
    path::AbstractString,
    grid_or_spectral_grid,
    dist_km;
    time_indices=Colon(),
    validate_initial_positions=true,
)
    """
    Compute FTLE from a SpeedyWeather `ParticleTracker` NetCDF file.

    Inputs:
        path: path to the particle-tracker NetCDF file
        grid_or_spectral_grid: spatial grid, or the SpectralGrid used for tracking
        dist_km: initial FTLE particle perturbation in km
        validate_initial_positions: validate the first particle positions
            against the canonical FTLE stencil before post-processing
        time_indices: optional time columns to process. Supports `:`, `:all`,
            `:first`, `:last`/`:final`, `:nonzero`/`:positive`, integer
            indices, integer-index iterables, and boolean masks.

    Outputs:
        FTLE_grid_time: FTLE matrix with one row per grid point and one column per recorded time
        time_hours: selected sampling times in hours

    Signed sampling times are treated as elapsed durations with magnitude
    `abs(time_hours)`. Zero-duration samples are returned as NaN because FTLE is
    undefined at T = 0. Non-finite sampling times throw `ArgumentError`.
    """

    npoints = _grid_npoints(grid_or_spectral_grid)
    particles_ds = NCDataset(path, "r")

    try
        _validate_particle_file(
            particles_ds,
            grid_or_spectral_grid,
            dist_km;
            validate_initial_positions,
        )
        time_hours = _particle_file_time_hours(particles_ds)
        _check_time_hours(time_hours)
        selected_time_indices = _checked_time_indices(time_indices, time_hours)
        FTLE_grid_time = Array{Float64}(undef, npoints, length(selected_time_indices))
        B = Array{Float64}(undef, 2, 2, npoints)

        FTLE_from_particles!(
            FTLE_grid_time,
            B,
            particles_ds["lon"],
            particles_ds["lat"],
            time_hours,
            npoints,
            dist_km;
            time_indices=selected_time_indices,
        )
        return FTLE_grid_time, _selected_time_hours(time_hours, selected_time_indices)
    finally
        close(particles_ds)
    end
end

export Re
export FTLE_from_particles!
export FTLE_from_particles
export FTLE_from_particle_file!
export FTLE_from_particle_file
