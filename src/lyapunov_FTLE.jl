using SpeedyWeather, RingGrids, NCDatasets

include("./FTLE_computations.jl")

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

function lyapunov_FTLE(; 
    simulation_days=10, 
    dist_km=10, 
    backwards=false, 
    dynamics=false, 
    model_type=BarotropicModel, 
    trunc::Int=31,
    spatial_grid_type=HEALPixGrid, 
    use_climatological=true,
    use_random=false,
    zonal_velocity_field::Field=zeros(HEALPixGrid(2)),
    meridional_velocity_field::Field=zeros(HEALPixGrid(2)),
    rint_hours=3
    )
    """
    TODO docstring is temporary
    Runs simulation with particle tracking and calculates the Finite-Time Lyapunov Exponent (FTLE).
    
    Inputs:
        simulation_days: number of days to run the simulation for
        dist_km: initial perturbation in km
        backwards: if true, run simulation backwards in time
        dynamics: if false, disables SpeedyWeather dynamics (makes velocity field static)
        model_type: model type. Currently only one-layer models are supported
        trunc: spectral truncation for the spectral grid
        spatial_grid_type: spatial grid type to use (e.g., HEALPixGrid, FullGaussianGrid)
        use_climatological: if true, use a fixed climatological velocity field
        use_random: if true, use a random velocity field
        zonal_velocity_field: Field representing the zonal velocity field (used if neither use_random nor use_climatological is true)
        meridional_velocity_field: Field representing the meridional velocity field (used if neither use_random nor use_climatological is true)
        rint_hours: sampling time for recording the positions of particles

    Outputs:
        FTLE_grid_time: NxM Matrix{Float64} . FTLE in units of 1/hour. N is number of grid points (spatial positions), M is number of time samples
        grid: grid object, that can tell you how indices in the first dimension of FTLE_grid_time map to a position on the sphere
        time_hours: Mx1 Vector{Float64}. Sampling times in hours
    """
    if use_random && use_climatological
        error("Cannot use both random and climatological velocity fields")
    end

    ### Setup simulation ###
    # Setup a default grid with 4 particles per grid cell
    spectral_grid = SpectralGrid(nlayers=1, trunc=trunc, Grid=spatial_grid_type)
    nparticles = 4*spectral_grid.npoints
    spectral_grid = SpectralGrid(nlayers=1, trunc=trunc, nparticles=nparticles, Grid=spatial_grid_type)

    # Set up particle advection scheme, model, and simulation
    particle_advection = ParticleAdvection2D(spectral_grid, backwards=backwards)
    model = model_type(spectral_grid, dynamics=dynamics; particle_advection)
    simulation = initialize!(model)

    # Apply prescribed velocity field
    if use_climatological
        # Set fixed velocity field to initial conditions of simulation
        # TODO don't understand the four lines of code below (they are from Milan)
        progn, diagn, model = SpeedyWeather.unpack(simulation)
        SpeedyWeather.scale!(progn, diagn, model.planet.radius)
        lf = 1
        transform!(diagn, progn, lf, model)
    elseif use_random
        u = randn(spectral_grid.grid)
        v = randn(spectral_grid.grid)
        simulation.diagnostic_variables.grid.u_grid .= u
        simulation.diagnostic_variables.grid.v_grid .= v
    else
        u_field = interpolate(spectral_grid.grid, zonal_velocity_field)
        v_field = interpolate(spectral_grid.grid, meridional_velocity_field)
        simulation.diagnostic_variables.grid.u_grid .= u_field
        simulation.diagnostic_variables.grid.v_grid .= v_field
    end

    # Add particle tracker to the model
    particle_tracker = ParticleTracker(spectral_grid, schedule=Schedule(every=Hour(rint_hours)))
    add!(model, :particle_tracker => particle_tracker)

    ### Perturb initial locations of particles ###
    londs, latds = RingGrids.get_londlatds(spectral_grid.grid)
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

    return FTLE_grid_time, spectral_grid.grid, time_hours

end

export lyapunov_FTLE
export Re