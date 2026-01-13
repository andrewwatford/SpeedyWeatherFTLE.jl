using SpeedyWeather, LinearAlgebra, GLMakie, NCDatasets

Re = 6.371e6 # Average Earth radius in meters

function setup_grid_FTLE(ntrunc, nlayers, grid_type)
    """
    Sets up a SpectralGrid with 4 particles per grid cell for calculating the FTLE

    ntrunc: truncation degree of the spectral grid
    nlayers: number of layers in the model atmosphere
    grid_type: grid type to use e.g. HEALPixGrid
    """
    spectral_grid = SpectralGrid(trunc=ntrunc, nlayers=nlayers, Grid=grid_type)
    nparticles = 4*spectral_grid.npoints # Number of particles
    spectral_grid = SpectralGrid(trunc=ntrunc, nlayers=nlayers, nparticles=nparticles, Grid=grid_type)
    return spectral_grid
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

function displacement_gradient_matrix_central(plonds, platds, dist_km)
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

    Ngpoints = length(plonds) ÷ 4 # Number of grid points

    dfac = Re / dist_km / 2

    cos_factor = cos.(deg2rad.(platds)) # Cos(latitude)

    # Derivative of x w.r.t. X
    xX = [deg2rad((plonds[4i-3] - plonds[4i-2]))*cos_factor[i]*dfac for i in 1:Ngpoints];
    # Derivative of y w.r.t. X 
    yX = [deg2rad((platds[4i-3] - platds[4i-2]))*dfac for i in 1:Ngpoints];
    # Derivative of x w.r.t. Y
    xY = [deg2rad((plonds[4i-1] - plonds[4i]))*cos_factor[i]*dfac for i in 1:Ngpoints];
    # Derivative of y w.r.t. Y 
    yY = [deg2rad((platds[4i-1] - platds[4i]))*dfac for i in 1:Ngpoints];

    # Displacement gradient matrix
    B = Array{eltype(xX)}(undef, 2, 2, Ngpoints)
    B[1,1,:] .= xX
    B[1,2,:] .= xY
    B[2,1,:] .= yX
    B[2,2,:] .= yY 
    # ^ Not the most elegant, but explicit 
    
    return B
end

function FTLE_from_eigenvalue(lmax, T)
    """
    Calculates the FTLE from a time interval and a maximum eigenvalue

    Inputs:
        lmax: maximum eigenvalue of the right Cauchy-Green deformation tensor
        T: time interval over which lmax was calculated. Units are reciprocal to units of T
    """
    return log(lmax)/2/T 
end

function FTLE_over_grid(B, T)
    """
    Compute FTLE over a grid

    Inputs:
        B: displacement gradient matrix
        T: time after particle release which B corresponds to

    Outputs: 
        FTLE_grid: FTLE at each grid point. Units are 
    """

    Ngpoints = size(B, 3) # Number of grid points

    FTLE_grid = Vector{Float64}(undef, Ngpoints) 

    for k in 1:Ngpoints
        # Gradient tensor
        Bk = B[:,:,k];
        # Right Cauchy-Green tensor
        CG = Bk' * Bk;
        # Largest eigenvalue
        lmax = maximum(eigvals(CG))
        # FTLE - in units of days^-1
        FTLE_grid[k] = FTLE_from_eigenvalue(lmax, T)
    end

    return FTLE_grid
end

function lyapunov_FTLE(; simulation_days=10, dist_km=10, backwards=false, 
    dynamics=false, grid_type=HEALPixGrid, model_type=BarotropicModel, use_initial=true,
    rint_hours=3)
    """
    TODO docstring is temporary
    Runs simulation with particle tracking and calculates the Finite-Time Lyapunov Exponent (FTLE)
    
    Inputs:
        simulation_days: number of days to run the simulation for
        dist_km: initial perturbation in km
        backwards: if true, run simulation backwards in time
        dynamics: if false, disables SpeedyWeather dynamics
        grid_type: spectral grid type
        model_type: model type. Currently only one-layer models are supported
        use_initial: if true, take initial velocity field of simulation. If false, take random initial velocity field
        rint_hours: sampling time for recording the positions of particles

    Outputs:
        FTLE_grid_time: NxM Matrix{Float64} . FTLE in units of 1/hour. N is number of grid points (spatial positions), M is number of time samples
        grid: grid object, that can tell you how indices in the first dimension of FTLE_grid_time map to a position on the sphere
        time_hours: Mx1 Vector{Float64}. Sampling times in hours
    """

    ### Setup simulation ###

    # TODO currently only supports 1-layer simulations

    # Setup spectral grid TODO as thing below is temporary
    ntrunc = 50 # Resolution (temporary)
    spectral_grid = setup_grid_FTLE(ntrunc, 1, grid_type)

    # Set up particle advection scheme, model, and simulation
    particle_advection = ParticleAdvection2D(spectral_grid, backwards=backwards)
    model = model_type(spectral_grid, dynamics=dynamics; particle_advection)
    simulation = initialize!(model)

    # Apply prescribed velocity field
    grid = spectral_grid.grid
    if use_initial # Set fixed velocity field to initial conditions of simulation
        # TODO don't understand the four lines of code below (they are from Milan)
        progn, diagn, model = SpeedyWeather.unpack(simulation)
        SpeedyWeather.scale!(progn, diagn, model.planet.radius)
        lf = 1
        transform!(diagn, progn, lf, model)
    else # Set fixed velocity field to random
        u = randn(grid)
        v = randn(grid)
        simulation.diagnostic_variables.grid.u_grid .= u
        simulation.diagnostic_variables.grid.v_grid .= v
    end
    # TODO add ability to prescribe velocity field

    # Add particle tracker
    particle_tracker = ParticleTracker(spectral_grid, schedule=Schedule(every=Hour(rint_hours)))
    add!(model, :particle_tracker => particle_tracker)


    ### Perturb initial locations of particles ###
    londs, latds = RingGrids.get_londlatds(grid)
    (; particles) = simulation.prognostic_variables
    perturb_positions_FTLE(particles, londs, latds, dist_km)


    ### Run ###
    run!(simulation, period=Day(simulation_days))


    ### Calculate time-dependent FTLE ###

    # Read in particle positions over time
    path = joinpath(model.output.run_folder, particle_tracker.filename) # Path to netCDF file
    particles_ds = NCDataset(path,"r")
    # TODO Would be nice to get output directly rather than write to then read from file. Cannot find how to do this in doc 

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

    return FTLE_grid_time, grid, time_hours

end
