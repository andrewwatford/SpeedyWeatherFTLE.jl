function _check_initial_FTLE_positions(londs, latds, dist_km)
    Npoints = length(londs)
    length(latds) == Npoints || throw(DimensionMismatch("londs and latds must have the same length"))
    dist_km > 0 || throw(ArgumentError("dist_km must be positive"))
    return Npoints, rad2deg(dist_km * 1000 / Re)
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
SpeedyWeatherFTLE's displacement-gradient reconstruction.
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
- `dynamics = false`: keep the prescribed velocity field static when `false`.
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
- `return_result = false`: return an [`FTLEResult`](@ref) instead of a tuple.
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
        dynamics: if false, disables SpeedyWeather dynamics (makes velocity field static)
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
    dist_km > 0 || throw(ArgumentError("dist_km must be positive"))
    particle_advection_every_n_time_steps >= 1 || throw(ArgumentError("particle_advection_every_n_time_steps must be at least 1"))
    rint_hours > 0 || throw(ArgumentError("rint_hours must be positive"))
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
    particle_advection = ParticleAdvection2D(
        spectral_grid;
        nparticles=n_particles,
        backwards=backwards,
        every_n_time_steps=particle_advection_every_n_time_steps,
    )
    if model_type != BarotropicModel
        @warn "get_FTLE currently only tested with BarotropicModel. Unexpected behaviour may occur."
    end
    model = with_logger(ConsoleLogger(stderr, Logging.Warn)) do
        model_type(spectral_grid; dynamics=dynamics, particle_advection=particle_advection)
    end
    simulation = with_logger(ConsoleLogger(stderr, Logging.Warn)) do
        initialize!(model)
    end

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
    with_logger(ConsoleLogger(stderr, Logging.Warn)) do
        SpeedyWeather.initialize!(simulation; period=Day(simulation_days), output=false)
        simulation.variables.grid.u[:, 1, 1] .= u
        simulation.variables.grid.v[:, 1, 1] .= v
        SpeedyWeather.initialize!(simulation.variables, particles, model)
        SpeedyWeather.time_stepping!(simulation)
        SpeedyWeather.finalize!(simulation)
    end

    ### Calculate time-dependent FTLE ###

    # Read in particle positions over time
    path = joinpath(particle_tracker.path == "" ? model.output.run_path : particle_tracker.path, particle_tracker.filename)
    should_keep_particle_file = keep_particle_file || return_particle_file_path

    try
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
export initial_FTLE_particle_positions!
export initial_FTLE_particle_positions
export positive_FTLE
export negative_FTLE
export Re
