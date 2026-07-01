@inline _wrapped_lon_diff(lond1, lond2) = mod(lond1 - lond2 + 180, 360) - 180

function displacement_gradient_matrix_central!(B, plonds, platds, dist_km)
    """
    Compute the displacement gradient matrix given particle positions at a fixed time
    Uses a central difference scheme to do this, with four points per grid cell

    Inputs:
        plonds: longitudes of particles
        platds: latitudes of particles
        dist_km: perturbation in position applied before starting simulation

    Outputs:
        B: (2,2,N) array where N is the number of grid points in the simulation.
        B[:,:,k] is the displacement gradient matrix for the k'th grid point. 
    """

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

@inline function _cauchy_green_invariants(a, b, c, d)
    C11 = a*a + c*c
    C12 = a*b + c*d
    C22 = b*b + d*d
    trC = C11 + C22
    discr = sqrt((C11 - C22)^2 + 4C12^2)
    return trC, discr
end

@inline function _cauchy_green_eigenvalues(a, b, c, d)
    trC, discr = _cauchy_green_invariants(a, b, c, d)
    λmax = (trC + discr) / 2
    λmin = max((trC - discr) / 2, zero(λmax))
    return λmin, λmax
end

@inline _largest_cauchy_green_eigenvalue(a, b, c, d) = last(_cauchy_green_eigenvalues(a, b, c, d))

"""
    cauchy_green_eigenvalues_over_grid!(λmin_grid, λmax_grid, B)

Compute the eigenvalues of the right Cauchy-Green deformation tensor `B'B` at
each grid point and write them to `λmin_grid` and `λmax_grid`.

This low-allocation diagnostic complements FTLE: `sqrt(λmax)` is the largest
finite-time stretching factor, while `sqrt(λmin)` is the smallest. Values are
clamped at zero only for the smaller eigenvalue to avoid tiny negative roundoff
from nearly singular deformations.
"""
function cauchy_green_eigenvalues_over_grid!(λmin_grid, λmax_grid, B)
    Ngpoints = size(B, 3)
    length(λmin_grid) == Ngpoints || throw(DimensionMismatch("λmin_grid must have length $Ngpoints"))
    length(λmax_grid) == Ngpoints || throw(DimensionMismatch("λmax_grid must have length $Ngpoints"))

    @inbounds for k in 1:Ngpoints
        λmin, λmax = _cauchy_green_eigenvalues(B[1, 1, k], B[1, 2, k], B[2, 1, k], B[2, 2, k])
        λmin_grid[k] = λmin
        λmax_grid[k] = λmax
    end

    return λmin_grid, λmax_grid
end

"""
    cauchy_green_eigenvalues_over_grid(B)

Allocating companion to [`cauchy_green_eigenvalues_over_grid!`](@ref). Returns
`λmin_grid, λmax_grid`.
"""
function cauchy_green_eigenvalues_over_grid(B)
    Ngpoints = size(B, 3)
    λmin_grid = Vector{Float64}(undef, Ngpoints)
    λmax_grid = Vector{Float64}(undef, Ngpoints)
    return cauchy_green_eigenvalues_over_grid!(λmin_grid, λmax_grid, B)
end

"""
    maximum_stretching_over_grid!(stretching_grid, B)

Compute the largest finite-time stretching factor `sqrt(λmax(B'B))` at each
grid point. This is equivalent to `exp(abs(T) * FTLE)` for nonzero integration
time `T`, but remains a purely kinematic deformation diagnostic that does not
require a time scale.
"""
function maximum_stretching_over_grid!(stretching_grid, B)
    Ngpoints = size(B, 3)
    length(stretching_grid) == Ngpoints || throw(DimensionMismatch("stretching_grid must have length $Ngpoints"))

    @inbounds for k in 1:Ngpoints
        λmax = _largest_cauchy_green_eigenvalue(B[1, 1, k], B[1, 2, k], B[2, 1, k], B[2, 2, k])
        stretching_grid[k] = sqrt(λmax)
    end

    return stretching_grid
end

"""
    maximum_stretching_over_grid(B)

Allocating companion to [`maximum_stretching_over_grid!`](@ref).
"""
function maximum_stretching_over_grid(B)
    stretching_grid = Vector{Float64}(undef, size(B, 3))
    return maximum_stretching_over_grid!(stretching_grid, B)
end

function FTLE_over_grid!(FTLE_grid, B, T)
    """
    Compute FTLE over a grid

    Inputs:
        B: displacement gradient matrix
        T: time after particle release which B corresponds to

    Outputs: 
        FTLE_grid: FTLE at each grid point. Units are 1 / [T].
        If T is zero, the FTLE is undefined and NaN is written.
    """

    Ngpoints = size(B, 3) # Number of grid points
    length(FTLE_grid) == Ngpoints || throw(DimensionMismatch("FTLE_grid must have length $Ngpoints"))

    if iszero(T)
        fill!(FTLE_grid, NaN)
        return FTLE_grid
    end

    twoT = 2 * T
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
    elseif time_indices in (:nonzero, :positive)
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
selected columns are determined by `time_indices`. `B` is a reusable work array
with dimensions `(2, 2, grid point)`.

Supported `time_indices` values are `:`, `:all`, `:first`, `:last`, `:final`,
`:nonzero`, `:positive`, an integer index, integer-index iterables, or a boolean
mask. Zero-duration samples are filled with `NaN` because FTLE is undefined at
`t = 0`.

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
        B: displacement-gradient work array with dimensions (2, 2, grid point)
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

    dist_km > 0 || throw(ArgumentError("dist_km must be positive"))
    size(plonds_time) == size(platds_time) ||
        throw(DimensionMismatch("plonds_time and platds_time must have the same size"))

    npoints = _grid_npoints(grid_or_npoints)
    npoints > 0 || throw(ArgumentError("number of grid points must be positive"))

    expected_nparticles = 4 * npoints
    actual_nparticles = size(plonds_time, 1)
    actual_nparticles == expected_nparticles ||
        throw(DimensionMismatch("particle trajectories have $actual_nparticles particles, expected $expected_nparticles"))

    n_times = length(time_hours)
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

    Zero-duration samples are returned as NaN because FTLE is undefined at T = 0.
    """

    npoints = _grid_npoints(grid_or_npoints)
    n_times = length(time_hours)
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
    time_vec = particles_ds["time"][:] # DateTime
    return Float64.((time_vec .- time_vec[1]) ./ Hour(1)) # Hours since release time
end

function _validate_particle_file(particles_ds, npoints)
    expected_nparticles = 4 * npoints
    actual_nparticles = particles_ds.dim["particle"]
    actual_nparticles == expected_nparticles ||
        throw(DimensionMismatch("particle file has $actual_nparticles particles, expected $expected_nparticles"))

    return nothing
end

"""
    FTLE_from_particle_file!(
        FTLE_grid_time,
        B,
        path,
        grid_or_spectral_grid,
        dist_km;
        time_indices = :
    )

Compute FTLE from a SpeedyWeather `ParticleTracker` NetCDF file, writing into
caller-provided arrays.

`path` must point to a particle-tracker file with longitude and latitude
variables named `lon` and `lat`. The file must contain four particles per FTLE
grid point. `grid_or_spectral_grid` may be the spatial grid or the `SpectralGrid`
used for the tracking run.

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
)
    """
    Compute FTLE from a SpeedyWeather `ParticleTracker` NetCDF file, reusing
    caller-provided output and work arrays.

    Inputs:
        FTLE_grid_time: output matrix with dimensions (grid point, time)
        B: displacement-gradient work array with dimensions (2, 2, grid point)
        path: path to the particle-tracker NetCDF file
        grid_or_spectral_grid: spatial grid, or the SpectralGrid used for tracking
        dist_km: initial FTLE particle perturbation in km
        time_indices: optional time columns to process. Supports `:`, `:all`,
            `:first`, `:last`/`:final`, `:nonzero`/`:positive`, integer
            indices, integer-index iterables, and boolean masks.

    Outputs:
        FTLE_grid_time: modified in place
        time_hours: selected sampling times in hours

    Zero-duration samples are returned as NaN because FTLE is undefined at T = 0.
    """

    npoints = _grid_npoints(grid_or_spectral_grid)
    particles_ds = NCDataset(path, "r")

    try
        _validate_particle_file(particles_ds, npoints)
        time_hours = _particle_file_time_hours(particles_ds)
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
    FTLE_from_particle_file(path, grid_or_spectral_grid, dist_km; time_indices = :)

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
)
    """
    Compute FTLE from a SpeedyWeather `ParticleTracker` NetCDF file.

    Inputs:
        path: path to the particle-tracker NetCDF file
        grid_or_spectral_grid: spatial grid, or the SpectralGrid used for tracking
        dist_km: initial FTLE particle perturbation in km
        time_indices: optional time columns to process. Supports `:`, `:all`,
            `:first`, `:last`/`:final`, `:nonzero`/`:positive`, integer
            indices, integer-index iterables, and boolean masks.

    Outputs:
        FTLE_grid_time: FTLE matrix with one row per grid point and one column per recorded time
        time_hours: selected sampling times in hours

    Zero-duration samples are returned as NaN because FTLE is undefined at T = 0.
    """

    npoints = _grid_npoints(grid_or_spectral_grid)
    particles_ds = NCDataset(path, "r")

    try
        _validate_particle_file(particles_ds, npoints)
        time_hours = _particle_file_time_hours(particles_ds)
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
export cauchy_green_eigenvalues_over_grid!
export cauchy_green_eigenvalues_over_grid
export maximum_stretching_over_grid!
export maximum_stretching_over_grid
export FTLE_from_particles!
export FTLE_from_particles
export FTLE_from_particle_file!
export FTLE_from_particle_file
