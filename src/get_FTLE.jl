function perturb_positions_FTLE(particles, londs, latds, dist_km)
    """
    Sets up the initial positions of particles for calculating the FTLE
    
    !! Modifies the simulation object in place

    particles: vector containing all Particle objects in the simulation
    londs: longitudes of grid cells
    latds: latitudes of grid cells
    dist_km: perturbation to apply in km
    """

    Npoints = length(londs) # Number of grid points

    cos_factor = cos.(deg2rad.(latds)) # Cos(latitude)

    del_lat = rad2deg((dist_km * 1000 / Re)) # Latitude perturbation in degrees

    # Perturbed East/West
    particles[1:4:end] .= [Particle(londs[i] + del_lat/cos_factor[i], latds[i]) for i in 1:Npoints]
    particles[2:4:end] .= [Particle(londs[i] - del_lat/cos_factor[i], latds[i]) for i in 1:Npoints]

    # Perturbed North/South
    particles[3:4:end] .= [Particle(londs[i], latds[i] + del_lat) for i in 1:Npoints]
    particles[4:4:end] .= [Particle(londs[i], latds[i] - del_lat) for i in 1:Npoints]
end

function get_FTLE(
    u::Field, 
    v::Field;
    simulation_days=10,
    dist_km=10,
    backwards=false,
    dynamics=false,
    rint_hours=3,
    model_type=BarotropicModel
)
    """
    TODO docstring is temporary
    Calculates the Finite-Time Lyapunov Exponent (FTLE) for a given velocity field.
    
    Inputs:
        u: Field representing the zonal velocity field
        v: Field representing the meridional velocity field
        simulation_days: number of days to run the simulation for
        dist_km: initial perturbation for released particles in km
        backwards: if true, run simulation backwards in time
        dynamics: if false, disables SpeedyWeather dynamics (makes velocity field static)
        rint_hours: sampling time for recording the positions of particles

    Outputs:
        FTLE_grid_time: NxM Matrix{Float64} . FTLE in units of 1/hour. N is number of grid points (spatial positions), M is number of time samples
        grid: grid object, that can tell you how indices in the first dimension of FTLE_grid_time map to a position on the sphere
        time_hours: Mx1 Vector{Float64}. Sampling times in hours
    """

    # Setup the spectral grid based on the input velocity fields
    if u.grid != v.grid
        error("Velocity fields u and v must be defined on the same grid")
    end
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
    n_particles = 4*temp_spectral_grid.npoints
    spectral_grid = SpectralGrid(nlayers=1, trunc=trunc, nparticles=n_particles, Grid=spatial_grid_type)

    # Set up particle advection scheme, model, and simulation
    particle_advection = ParticleAdvection2D(spectral_grid, backwards=backwards)
    if model_type != BarotropicModel
        @warn "get_FTLE currently only tested with BarotropicModel. Unexpected behaviour may occur."
    end
    model = model_type(spectral_grid, dynamics=dynamics; particle_advection=particle_advection)
    simulation = initialize!(model)

    # Apply prescribed velocity field
    simulation.diagnostic_variables.grid.u_grid .= u
    simulation.diagnostic_variables.grid.v_grid .= v

    # Add particle tracker to the model
    particle_tracker = ParticleTracker(spectral_grid, schedule=Schedule(every=Hour(rint_hours)))
    add!(model, :particle_tracker => particle_tracker)

    ### Perturb initial locations of particles ###
    londs, latds = RingGrids.get_londlatds(spatial_grid)
    (; particles) = simulation.prognostic_variables
    perturb_positions_FTLE(particles, londs, latds, dist_km)

    ### Run ###
    run!(simulation, period=Day(simulation_days))

    ### Calculate time-dependent FTLE ###

    # Read in particle positions over time
    path = joinpath(model.output.run_folder, particle_tracker.filename) # Path to netCDF file
    particles_ds = NCDataset(path,"r")

    # Time dimension
    time_vec = particles_ds["time"][:] # DateTime
    time_hours = ( time_vec .- time_vec[1] ) ./ Hour(1) # Hours since release time

    # Initialise array to hold FTLE over grid and over time
    FTLE_grid_time = Array{Float64}(undef, spectral_grid.npoints, particles_ds.dim["time"]) 

    # Iterate over time steps
    for (tindex, thour) in enumerate(time_hours)

        # Particle latitudes and longitudes
        plonds = particles_ds["lon"][:,tindex]
        platds = particles_ds["lat"][:,tindex]

        # Calculate displacement gradient matrix
        B = displacement_gradient_matrix_central(plonds, platds, dist_km)

        # Calculate FTLE
        FTLE_grid_time[:,tindex] .= FTLE_over_grid(B, thour)

    end

    close(particles_ds)
    rm(path) # Remove temporary netCDF file

    return FTLE_grid_time, spectral_grid, time_hours

end

export get_FTLE
export Re