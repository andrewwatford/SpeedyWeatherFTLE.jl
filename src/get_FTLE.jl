const _MAX_LOCAL_STENCIL_DEGREES = 20.0

const _SPEEDYWEATHER_SOURCE_INSTALL_MESSAGE =
    "SpeedyWeatherFTLE currently requires the mk/lyapunov2 branch of the " *
    "SpeedyWeather monorepo. Registered SpeedyWeather releases do not yet " *
    "provide the ParticleAdvection2D keyword API required by get_FTLE. " *
    "Install the SpeedyWeather subpackages from " *
    "https://github.com/SpeedyWeather/SpeedyWeather.jl at rev mk/lyapunov2 " *
    "before adding SpeedyWeatherFTLE."

function _ftle_particle_advection_2d(spectral_grid; nparticles, backwards, every_n_time_steps)
    try
        return ParticleAdvection2D(
            spectral_grid;
            nparticles,
            backwards,
            every_n_time_steps,
        )
    catch err
        if err isa MethodError
            throw(ArgumentError(_SPEEDYWEATHER_SOURCE_INSTALL_MESSAGE))
        end
        rethrow()
    end
end

function _check_initial_FTLE_positions(londs, latds, dist_km)
    Npoints = length(londs)
    length(latds) == Npoints || throw(DimensionMismatch("londs and latds must have the same length"))
    _check_dist_km(dist_km)

    del_lat = rad2deg(dist_km * 1000 / Re)
    del_lat <= _MAX_LOCAL_STENCIL_DEGREES ||
        throw(ArgumentError("dist_km=$dist_km creates a north/south FTLE stencil offset of $(del_lat)°, which is too large for the local spherical finite-difference approximation"))

    @inbounds for i in eachindex(londs, latds)
        isfinite(londs[i]) && isfinite(latds[i]) ||
            throw(ArgumentError("FTLE stencil center positions must be finite"))
        abs(latds[i]) < 90 ||
            throw(ArgumentError("FTLE stencil center latitude $(latds[i])° at index $i is at or beyond a pole"))
        abs(latds[i]) + del_lat < 90 ||
            throw(ArgumentError("dist_km=$dist_km makes the FTLE stencil cross a pole at latitude $(latds[i])° (index $i)"))

        del_lon = del_lat / abs(cosd(latds[i]))
        isfinite(del_lon) && del_lon <= _MAX_LOCAL_STENCIL_DEGREES ||
            throw(ArgumentError("dist_km=$dist_km creates an east/west FTLE stencil offset of $(del_lon)° at latitude $(latds[i])° (index $i), which is too large for the local spherical finite-difference approximation"))
    end

    return Npoints, del_lat
end

"""
    initial_FTLE_particle_positions!(plonds, platds, londs, latds, dist_km)
    initial_FTLE_particle_positions!(plonds, platds, grid_or_spectral_grid, dist_km)

Fill longitude and latitude arrays with the canonical four-particle FTLE
release stencil around each grid point.

Particles are written in east, west, north, south order for each grid point,
matching [`FTLE_from_particles`](@ref), [`FTLE_from_particle_file`](@ref), and
SpeedyWeather's `ParticleTracker` post-processing. `plonds` and `platds` must
each contain four entries per grid point.
"""
function initial_FTLE_particle_positions!(plonds, platds, londs, latds, dist_km)
    Npoints, del_lat = _check_initial_FTLE_positions(londs, latds, dist_km)
    length(plonds) == 4 * Npoints || throw(DimensionMismatch("plonds must contain four entries per grid point"))
    length(platds) == 4 * Npoints || throw(DimensionMismatch("platds must contain four entries per grid point"))

    @inbounds for i in 1:Npoints
        p = 4i
        del_lon = del_lat / cosd(latds[i])

        plonds[p - 3] = londs[i] + del_lon
        platds[p - 3] = latds[i]

        plonds[p - 2] = londs[i] - del_lon
        platds[p - 2] = latds[i]

        plonds[p - 1] = londs[i]
        platds[p - 1] = latds[i] + del_lat

        plonds[p] = londs[i]
        platds[p] = latds[i] - del_lat
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

Return longitude and latitude vectors for the canonical four-particle FTLE
release stencil around each grid point.

The returned vectors have length `4length(londs)`. For each grid point, entries
are ordered east, west, north, south, which is the layout expected by
[`FTLE_from_particles`](@ref), [`FTLE_from_particle_file`](@ref), and
SpeedyWeatherFTLE's flow-map Jacobian reconstruction. The stencil is a local
spherical finite-difference approximation; requests with non-finite or
non-positive `dist_km`, pole-crossing particles, high-latitude singular
east/west offsets, or clearly nonlocal separations throw `ArgumentError`.
"""
function initial_FTLE_particle_positions(londs, latds, dist_km)
    Npoints = length(londs)
    plonds = Vector{Float64}(undef, 4 * Npoints)
    platds = Vector{Float64}(undef, 4 * Npoints)
    return initial_FTLE_particle_positions!(plonds, platds, londs, latds, dist_km)
end

function initial_FTLE_particle_positions(grid_or_spectral_grid, dist_km)
    npoints = _grid_npoints(grid_or_spectral_grid)
    plonds = Vector{Float64}(undef, 4 * npoints)
    platds = Vector{Float64}(undef, 4 * npoints)
    return initial_FTLE_particle_positions!(plonds, platds, grid_or_spectral_grid, dist_km)
end

function perturb_positions_FTLE(particles, londs, latds, dist_km)
    """
    Sets up the initial positions of particles for calculating the FTLE
    
    !! Modifies the simulation object in place

    particles: vector containing all Particle objects in the simulation
    londs: longitudes of grid cells
    latds: latitudes of grid cells
    dist_km: perturbation to apply in km
    """

    Npoints, del_lat = _check_initial_FTLE_positions(londs, latds, dist_km)
    length(particles) == 4 * Npoints || throw(DimensionMismatch("particles must contain four particles per grid point"))

    @inbounds for i in 1:Npoints
        p = 4i
        del_lon = del_lat / cosd(latds[i])

        particles[p - 3] = Particle(londs[i] + del_lon, latds[i])
        particles[p - 2] = Particle(londs[i] - del_lon, latds[i])
        particles[p - 1] = Particle(londs[i], latds[i] + del_lat)
        particles[p] = Particle(londs[i], latds[i] - del_lat)
    end
end

const _FTLE_PARTICLE_ORDER_SYMBOLS = (:east, :west, :north, :south)

"""
    FTLEParticleSetup

Metadata returned by [`prepare_FTLE_particles!`](@ref) and
[`attach_FTLE_tracker!`](@ref).

The setup records the grid and `dist_km` needed for safe post-processing, the
particle order expected by the finite-difference reconstruction, the saved
particle file path when a tracker has been attached, and the default
`time_indices` selector to use when converting the saved trajectory file into
FTLE values.
"""
struct FTLEParticleSetup{S, G, T, P, I}
    spectral_grid::S
    grid::G
    dist_km::Float64
    particle_order::NTuple{4, Symbol}
    particle_file_path::P
    tracker::T
    callback_name::Union{Nothing, Symbol}
    particle_tracker_keepbits::Union{Nothing, Int}
    time_semantics::Symbol
    time_indices::I
end

function _simulation_model(simulation)
    hasproperty(simulation, :model) ||
        throw(ArgumentError("expected a SpeedyWeather Simulation with a model field"))
    return simulation.model
end

function _simulation_spectral_grid(simulation)
    model = _simulation_model(simulation)
    hasproperty(model, :spectral_grid) ||
        throw(ArgumentError("simulation.model must expose a spectral_grid"))
    return model.spectral_grid
end

function _ftle_setup(
    spectral_grid,
    dist_km;
    particle_file_path=nothing,
    tracker=nothing,
    callback_name=nothing,
    particle_tracker_keepbits=nothing,
    time_indices=Colon(),
)
    return FTLEParticleSetup(
        spectral_grid,
        _spatial_grid(spectral_grid),
        Float64(dist_km),
        _FTLE_PARTICLE_ORDER_SYMBOLS,
        particle_file_path,
        tracker,
        callback_name,
        particle_tracker_keepbits,
        :elapsed_hours_since_release,
        time_indices,
    )
end

"""
    prepare_FTLE_particles!(simulation; dist_km = 10, spectral_grid = simulation.model.spectral_grid)

Mutate an existing SpeedyWeather simulation's particles into the canonical FTLE
release stencil.

The simulation must already have a particle-advection scheme with four
particles per grid point. The helper writes particles in east, west, north,
south order and returns an [`FTLEParticleSetup`](@ref) containing the grid,
`dist_km`, particle order, and time semantics needed for later
post-processing.
"""
function prepare_FTLE_particles!(
    simulation;
    dist_km=10,
    spectral_grid=_simulation_spectral_grid(simulation),
    time_indices=Colon(),
)
    _check_dist_km(dist_km)
    londs, latds = RingGrids.get_londlatds(_spatial_grid(spectral_grid))
    hasproperty(simulation, :variables) &&
        hasproperty(simulation.variables, :prognostic) &&
        hasproperty(simulation.variables.prognostic, :particles) ||
        throw(ArgumentError("simulation.variables.prognostic.particles is required to prepare FTLE particles"))

    perturb_positions_FTLE(simulation.variables.prognostic.particles, londs, latds, dist_km)
    return _ftle_setup(spectral_grid, dist_km; time_indices)
end

"""
    attach_FTLE_tracker!(simulation; kwargs...)

Prepare FTLE particles, attach a SpeedyWeather `ParticleTracker` callback, and
return metadata for post-processing.

The simulation's model must already include a particle-advection scheme with
four particles per grid point. By default, particle output is written to a fresh
temporary directory so the returned `particle_file_path` is stable. Pass an
explicit `path` when you want to keep the file in a known location.

# Keyword Arguments

- `dist_km = 10`: FTLE particle perturbation distance in kilometres.
- `rint_hours = 3`: particle output cadence in hours.
- `path = mktempdir()`: directory for the particle file; must be non-empty.
- `filename = "particles.nc"`: particle file name.
- `particle_tracker_keepbits = 15`: mantissa bits retained in NetCDF output.
- `particle_tracker_compression_level = 1`: NetCDF compression level, from
  `0` to `9`.
- `particle_tracker_shuffle = false`: enable the NetCDF shuffle filter.
- `callback_name = :ftle_particle_tracker`: callback key in
  `simulation.model.callbacks`.
- `replace = false`: pass `true` to replace an existing callback with the same
  name.
- `time_indices = :`: default time selector used by
  `FTLE_from_particle_file(setup)`.
- `prepare_particles = true`: set to `false` only if you have already called
  [`prepare_FTLE_particles!`](@ref) with matching metadata.
"""
function attach_FTLE_tracker!(
    simulation;
    dist_km=10,
    rint_hours=3,
    path=mktempdir(; prefix="speedyweatherftle-particles-"),
    filename="particles.nc",
    particle_tracker_keepbits=15,
    particle_tracker_compression_level=1,
    particle_tracker_shuffle=false,
    callback_name::Symbol=:ftle_particle_tracker,
    replace::Bool=false,
    time_indices=Colon(),
    prepare_particles::Bool=true,
)
    _check_dist_km(dist_km)
    _check_positive_hours(rint_hours, "rint_hours")
    particle_tracker_keepbits >= 1 || throw(ArgumentError("particle_tracker_keepbits must be positive"))
    0 <= particle_tracker_compression_level <= 9 ||
        throw(ArgumentError("particle_tracker_compression_level must be between 0 and 9"))
    isempty(path) &&
        throw(ArgumentError("attach_FTLE_tracker! requires a non-empty path so particle_file_path is known"))
    isempty(filename) && throw(ArgumentError("filename must be non-empty"))

    model = _simulation_model(simulation)
    spectral_grid = _simulation_spectral_grid(simulation)

    if !replace && haskey(model.callbacks, callback_name)
        throw(ArgumentError("simulation.model.callbacks already has callback $callback_name; pass replace=true to replace it"))
    end

    mkpath(path)
    if prepare_particles
        prepare_FTLE_particles!(simulation; dist_km, spectral_grid, time_indices)
    end

    particle_tracker = ParticleTracker(
        spectral_grid;
        schedule=Schedule(every=Hour(rint_hours)),
        keepbits=particle_tracker_keepbits,
        compression_level=particle_tracker_compression_level,
        shuffle=particle_tracker_shuffle,
        path,
        filename,
    )
    model.callbacks[callback_name] = particle_tracker

    return _ftle_setup(
        spectral_grid,
        dist_km;
        particle_file_path=joinpath(path, filename),
        tracker=particle_tracker,
        callback_name,
        particle_tracker_keepbits=Int(particle_tracker_keepbits),
        time_indices,
    )
end

function _particle_setup_path(setup::FTLEParticleSetup)
    setup.particle_file_path !== nothing ||
        throw(ArgumentError("FTLEParticleSetup has no particle_file_path; use attach_FTLE_tracker! before file post-processing"))
    isfile(setup.particle_file_path) ||
        throw(ArgumentError("particle file $(setup.particle_file_path) does not exist; run and finalize the SpeedyWeather simulation before post-processing"))
    return setup.particle_file_path
end

function _write_particle_setup_metadata(setup::FTLEParticleSetup)
    path = _particle_setup_path(setup)
    keepbits = setup.particle_tracker_keepbits === nothing ? 15 : setup.particle_tracker_keepbits
    _write_FTLE_particle_file_metadata(path; dist_km=setup.dist_km, particle_tracker_keepbits=keepbits)
    return path
end

"""
    FTLE_from_particle_file(setup::FTLEParticleSetup; time_indices = setup.time_indices, validate_initial_positions = true)

Post-process the particle file recorded by [`attach_FTLE_tracker!`](@ref).

This convenience method writes SpeedyWeatherFTLE compatibility metadata to the
saved particle file before calling the path-based
[`FTLE_from_particle_file`](@ref).
"""
function FTLE_from_particle_file(
    setup::FTLEParticleSetup;
    time_indices=setup.time_indices,
    validate_initial_positions=true,
)
    path = _write_particle_setup_metadata(setup)
    return FTLE_from_particle_file(
        path,
        setup.spectral_grid,
        setup.dist_km;
        time_indices,
        validate_initial_positions,
    )
end

"""
    FTLE_from_particle_file!(FTLE_grid_time, B, setup::FTLEParticleSetup; kwargs...)

In-place companion to [`FTLE_from_particle_file`](@ref) for an
[`FTLEParticleSetup`](@ref).
"""
function FTLE_from_particle_file!(
    FTLE_grid_time,
    B,
    setup::FTLEParticleSetup;
    time_indices=setup.time_indices,
    validate_initial_positions=true,
)
    path = _write_particle_setup_metadata(setup)
    return FTLE_from_particle_file!(
        FTLE_grid_time,
        B,
        path,
        setup.spectral_grid,
        setup.dist_km;
        time_indices,
        validate_initial_positions,
    )
end

"""
    get_FTLE(u::Field, v::Field; kwargs...)

Run a SpeedyWeather particle-tracking simulation from prescribed zonal and
meridional velocity fields, then compute finite-time Lyapunov exponents from
the tracked particle trajectories.

`u` and `v` must be `RingGrids.Field` objects on the same grid. Four particles
are released around each grid point with initial separation `dist_km`, the
trajectories are written by SpeedyWeather's `ParticleTracker`, and the saved
particle positions are post-processed with [`FTLE_from_particle_file`](@ref).

# Main Keyword Arguments

- `simulation_days = 10`: simulation duration in days.
- `dist_km = 10`: particle perturbation distance in kilometres.
- `backwards = false`: run backward in time for negative-time FTLE.
- `dynamics = false`: required. The prescribed-field wrapper supports frozen
  supplied velocity fields only; `dynamics = true` is rejected because it does
  not yet initialize SpeedyWeather prognostic state from `u`/`v`.
- `rint_hours = 3`: particle output cadence in hours.
- `model_type = BarotropicModel`: SpeedyWeather model type.
- `particle_advection_every_n_time_steps = 6`: particle advection cadence.
- `particle_tracker_keepbits = 15`: mantissa bits retained in NetCDF output.
- `particle_tracker_compression_level = 1`: NetCDF compression level, from `0`
  to `9`.
- `particle_tracker_shuffle = false`: enable the NetCDF shuffle filter.
- `particle_tracker_path = ""`: directory for the particle file.
- `particle_tracker_filename = "particles.nc"`: particle file name.
- `keep_particle_file = false`: keep the particle file after FTLE computation.
- `return_particle_file_path = false`: return the particle file path and keep
  the file.
- `return_result = false`: when true, return an [`FTLEResult`](@ref) instead
  of a tuple.
- `time_indices = :`: particle-tracker output columns to post-process.

Supported `time_indices` values are `:`, `:all`, `:first`, `:last`, `:final`,
`:nonzero`, `:positive`, an integer index, integer-index iterables, or a boolean
mask.

# Returns

By default, returns `FTLE_grid_time, spectral_grid, time_hours`.

When `return_particle_file_path = true`, returns
`FTLE_grid_time, spectral_grid, time_hours, particle_file_path` and keeps the
particle file.

When `return_result = true`, returns an [`FTLEResult`](@ref) with the same data
and run metadata.
"""
function get_FTLE(
    u::Field, 
    v::Field;
    simulation_days=10,
    dist_km=10,
    backwards=false,
    dynamics=false,
    rint_hours=3,
    model_type=BarotropicModel,
    particle_advection_every_n_time_steps=6,
    particle_tracker_keepbits=15,
    particle_tracker_compression_level=1,
    particle_tracker_shuffle=false,
    particle_tracker_path="",
    particle_tracker_filename="particles.nc",
    keep_particle_file=false,
    return_particle_file_path=false,
    return_result=false,
    time_indices=Colon(),
)
    """
    Calculates the Finite-Time Lyapunov Exponent (FTLE) for a given velocity field.
    
    Inputs:
        u: Field representing the zonal velocity field
        v: Field representing the meridional velocity field
        simulation_days: number of days to run the simulation for
        dist_km: initial perturbation for released particles in km
        backwards: if true, run simulation backwards in time
        dynamics: must be false; the prescribed-field wrapper supports frozen velocity fields only
        rint_hours: sampling time for recording the positions of particles
        particle_advection_every_n_time_steps: advect particles every n model timesteps
        particle_tracker_keepbits: mantissa bits retained when particle positions are written to netCDF
        particle_tracker_compression_level: netCDF compression level for particle trajectories
        particle_tracker_shuffle: whether to use the netCDF shuffle filter for particle trajectories
        particle_tracker_path: directory for the temporary particle-tracker NetCDF file
        particle_tracker_filename: file name for the particle-tracker NetCDF file
        keep_particle_file: if true, do not delete the particle-tracker NetCDF file after computing FTLE
        return_particle_file_path: if true, also return the particle-tracker NetCDF path and keep the file
        return_result: if true, return an FTLEResult with named fields and metadata
        time_indices: optional particle-tracker time columns to post-process.
            Supports `:`, `:all`, `:first`, `:last`/`:final`,
            `:nonzero`/`:positive`, integer indices, integer-index iterables,
            and boolean masks.

    Outputs:
        FTLE_grid_time: NxM Matrix{Float64} . FTLE in units of 1/hour. N is number of grid points (spatial positions), M is number of time samples
        grid: grid object, that can tell you how indices in the first dimension of FTLE_grid_time map to a position on the sphere
        time_hours: Mx1 Vector{Float64}. Selected sampling times in hours
        particle_file_path: returned as a fourth value only when return_particle_file_path is true
        FTLEResult: returned instead of the tuple when return_result is true
    """

    # Setup the spectral grid based on the input velocity fields
    if u.grid != v.grid
        error("Velocity fields u and v must be defined on the same grid")
    end
    _check_dist_km(dist_km)
    if dynamics
        throw(ArgumentError("dynamics=true is not supported by the prescribed-field get_FTLE(u, v) wrapper: it does not yet initialize SpeedyWeather prognostic state from the supplied u/v fields; use dynamics=false for frozen-flow particle tracking"))
    end
    particle_advection_every_n_time_steps >= 1 || throw(ArgumentError("particle_advection_every_n_time_steps must be at least 1"))
    _check_positive_hours(rint_hours, "rint_hours")
    particle_tracker_keepbits >= 1 || throw(ArgumentError("particle_tracker_keepbits must be positive"))
    0 <= particle_tracker_compression_level <= 9 ||
        throw(ArgumentError("particle_tracker_compression_level must be between 0 and 9"))

    spatial_grid = u.grid
    spatial_grid_type = typeof(spatial_grid)
    # Since we do not give the user the ability to control dealiasing or truncation,
        # we use the default dealiasing of 2 and set truncation accordingly
    J = length(spatial_grid.rings)
    # Guess truncation from number of latitudinal rings
    T_lower = convert(Int, floor(2 * J / 3))
    T_upper = convert(Int, ceil(2 * J / 3))
    sg_lower = SpectralGrid(nlayers=1, trunc=T_lower, Grid=spatial_grid_type)
    sg_upper = SpectralGrid(nlayers=1, trunc=T_upper, Grid=spatial_grid_type)
    if spatial_grid == sg_lower.grid
        trunc = T_lower
    elseif spatial_grid == sg_upper.grid
        trunc = T_upper
    else
        error("Could not determine spectral truncation from provided grid")
    end
    temp_spectral_grid = SpectralGrid(nlayers=1, trunc=trunc, Grid=spatial_grid_type)
    n_particles = 4 * temp_spectral_grid.npoints
    spectral_grid = SpectralGrid(nlayers=1, trunc=trunc, Grid=spatial_grid_type)

    # Set up particle advection scheme, model, and simulation
    particle_advection = _ftle_particle_advection_2d(
        spectral_grid;
        nparticles=n_particles,
        backwards=backwards,
        every_n_time_steps=particle_advection_every_n_time_steps,
    )
    if model_type != BarotropicModel
        @warn "get_FTLE currently only tested with BarotropicModel. Unexpected behaviour may occur."
    end
    model = model_type(spectral_grid; dynamics=dynamics, particle_advection=particle_advection)
    simulation = initialize!(model)

    particle_tracker = ParticleTracker(
        spectral_grid;
        schedule=Schedule(every=Hour(rint_hours)),
        keepbits=particle_tracker_keepbits,
        compression_level=particle_tracker_compression_level,
        shuffle=particle_tracker_shuffle,
        path=particle_tracker_path,
        filename=particle_tracker_filename,
    )
    model.callbacks[:particle_tracker] = particle_tracker

    ### Perturb initial locations of particles ###
    londs, latds = RingGrids.get_londlatds(spatial_grid)
    (; particles) = simulation.variables.prognostic
    perturb_positions_FTLE(particles, londs, latds, dist_km)

    # SpeedyWeather.initialize!(simulation) transforms model prognostics to grid space,
    # so apply the prescribed static grid velocities after that initialization step.
    SpeedyWeather.initialize!(simulation; period=Day(simulation_days))
    simulation.variables.grid.u[:, 1, 1] .= u
    simulation.variables.grid.v[:, 1, 1] .= v
    SpeedyWeather.initialize!(simulation.variables, particles, model)
    SpeedyWeather.time_stepping!(simulation)
    SpeedyWeather.finalize!(simulation)

    ### Calculate time-dependent FTLE ###

    # Read in particle positions over time
    path = joinpath(particle_tracker.path == "" ? model.output.run_path : particle_tracker.path, particle_tracker.filename)
    should_keep_particle_file = keep_particle_file || return_particle_file_path

    try
        _write_FTLE_particle_file_metadata(path; dist_km, particle_tracker_keepbits)
        FTLE_grid_time, time_hours = FTLE_from_particle_file(path, spectral_grid, dist_km; time_indices)
        if return_result
            return FTLEResult(
                FTLE_grid_time,
                spectral_grid,
                time_hours;
                particle_file_path=should_keep_particle_file ? path : nothing,
                dist_km,
                backwards,
                dynamics,
                rint_hours,
            )
        elseif return_particle_file_path
            return FTLE_grid_time, spectral_grid, time_hours, path
        else
            return FTLE_grid_time, spectral_grid, time_hours
        end
    finally
        should_keep_particle_file || rm(path; force=true) # Remove temporary netCDF file
    end

end

function _reject_backwards_keyword(kwargs, function_name)
    if :backwards in keys(kwargs)
        throw(ArgumentError("$function_name fixes the FTLE time direction; use get_FTLE to pass backwards explicitly"))
    end
    return nothing
end

"""
    positive_FTLE(u::Field, v::Field; kwargs...)

Compute positive-time FTLE by calling [`get_FTLE`](@ref) with
`backwards = false`.

All other keyword arguments are forwarded to [`get_FTLE`](@ref). Passing a
`backwards` keyword is rejected because this wrapper fixes the time direction.
"""
function positive_FTLE(u::Field, v::Field; kwargs...)
    """
    Compute positive-time FTLE by running `get_FTLE` with `backwards=false`.
    """
    _reject_backwards_keyword(kwargs, "positive_FTLE")
    return get_FTLE(u, v; backwards=false, kwargs...)
end

"""
    negative_FTLE(u::Field, v::Field; kwargs...)

Compute negative-time FTLE by calling [`get_FTLE`](@ref) with
`backwards = true`.

All other keyword arguments are forwarded to [`get_FTLE`](@ref). Passing a
`backwards` keyword is rejected because this wrapper fixes the time direction.
"""
function negative_FTLE(u::Field, v::Field; kwargs...)
    """
    Compute negative-time FTLE by running `get_FTLE` with `backwards=true`.
    """
    _reject_backwards_keyword(kwargs, "negative_FTLE")
    return get_FTLE(u, v; backwards=true, kwargs...)
end

export get_FTLE
export FTLEParticleSetup
export initial_FTLE_particle_positions!
export initial_FTLE_particle_positions
export prepare_FTLE_particles!
export attach_FTLE_tracker!
export positive_FTLE
export negative_FTLE
export Re
