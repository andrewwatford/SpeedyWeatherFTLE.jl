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
    FTLE(u::Field, v::Field; kwargs...)

Compute FTLE for frozen zonal and meridional velocity fields, then return
`(ftle, spectral_grid, time_hours)`.
"""
function FTLE(
    u::Field,
    v::Field;
    simulation_days=10,
    dist_km=10,
    backwards=false,
    dynamics=false,
    rint_hours=3,
    particle_advection_every_n_time_steps=6,
    particle_tracker_keepbits=15,
    return_result=false,
    time_indices=Colon(),
)
    u.grid == v.grid || throw(ArgumentError("u and v must use the same grid"))
    dist_km = _check_dist_km(dist_km)
    dynamics && throw(ArgumentError("FTLE supports frozen prescribed fields only; pass dynamics=false"))
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
    )

    londs, latds = RingGrids.get_londlatds(u.grid)
    _perturb_positions!(simulation.variables.prognostic.particles, londs, latds, dist_km)

    SpeedyWeather.initialize!(simulation; period=_period_days(simulation_days))
    simulation.variables.grid.u[:, 1, 1] .= u
    simulation.variables.grid.v[:, 1, 1] .= v
    SpeedyWeather.initialize!(simulation.variables, simulation.variables.prognostic.particles, model)
    SpeedyWeather.time_stepping!(simulation)
    SpeedyWeather.finalize!(simulation)

    tracker = model.callbacks[:particle_tracker]
    path = joinpath(tracker.path == "" ? model.output.run_path : tracker.path, tracker.filename)

    try
        _write_particle_file_metadata(path; dist_km)
        ftle, time_hours = _FTLE_from_particle_file(path, spectral_grid, dist_km; time_indices)
        if return_result
            return FTLEResult(
                ftle,
                spectral_grid,
                time_hours;
                dist_km,
                backwards,
                dynamics=false,
                rint_hours,
            )
        else
            return ftle, spectral_grid, time_hours
        end
    finally
        rm(path; force=true)
    end
end
